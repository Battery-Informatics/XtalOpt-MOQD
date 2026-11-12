/**********************************************************************
  MomeTest - Unit tests for the MOME-X archive math (search/mome.h)

  Copyright (C) 2026 Samad Hajinazar

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
 ***********************************************************************/

#include <search/mome.h>
#include <search/selection.h>

#include <QtTest>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <set>

using namespace Search;

namespace {

MomeTessellation grid2d(int xBins, int yBins)
{
  MomeTessellation t;
  t.kind = MomeTessellation::Grid;
  t.lo = { 0.0, 10.0 };
  t.hi = { 1.0, 20.0 };
  t.bins = { xBins, yBins };
  return t;
}

int cellOf(const MomeTessellation& t, double a, double b)
{
  const double d[2] = { a, b };
  return assignCell(d, t);
}

} // namespace

class MomeTest : public QObject
{
  Q_OBJECT

private slots:
  // Tessellation
  void gridIndexIsRowMajor();
  void gridClampsOutOfRangeValues();
  void gridRejectsNonFiniteAndInvalidTessellation();
  void gridCoordsRoundTrip();
  void gridSupportsOneAndThreeDimensions();
  void cvtCentroidsAreDeterministicForSeed();
  void cvtCentroidsSpreadOverUnitBox();
  void cvtAssignsNearestCentroid();

  // Per-cell fronts
  void truncationKeepsBoundaryAndDropsMostCrowded();
  void truncationLeavesSmallFrontsAlone();

  // Sampling
  void weightsReplaceInfiniteWithMaxFinite();
  void weightsUniformForSmallOrDegenerateFronts();
  void samplingFrequenciesFollowCrowding();

  // Archive
  void archiveKeepsFirstFrontPerCell();
  void archiveRespectsMaxFrontSize();
  void archiveSkipsUnassignableRows();
};

void MomeTest::gridIndexIsRowMajor()
{
  const MomeTessellation t = grid2d(4, 5);
  QCOMPARE(t.cellCount(), 20LL);
  QCOMPARE(cellOf(t, 0.0, 10.0), 0);
  QCOMPARE(cellOf(t, 0.0, 12.5), 1);   // second bin of the last dimension
  QCOMPARE(cellOf(t, 0.30, 10.0), 5);  // second bin of the first dimension
  QCOMPARE(cellOf(t, 0.99, 19.9), 19);
  QCOMPARE(cellOf(t, 0.25, 12.0), 6);  // bin edges belong to the upper bin
}

void MomeTest::gridClampsOutOfRangeValues()
{
  const MomeTessellation t = grid2d(4, 5);
  QCOMPARE(cellOf(t, 1.0, 20.0), 19);   // upper bound goes to the last bin
  QCOMPARE(cellOf(t, 7.0, 99.0), 19);
  QCOMPARE(cellOf(t, -3.0, -3.0), 0);
  QCOMPARE(cellOf(t, -3.0, 99.0), 4);
}

void MomeTest::gridRejectsNonFiniteAndInvalidTessellation()
{
  MomeTessellation t = grid2d(4, 5);
  QCOMPARE(cellOf(t, std::numeric_limits<double>::quiet_NaN(), 12.0), -1);
  QCOMPARE(cellOf(t, 0.5, std::numeric_limits<double>::infinity()), -1);

  t.hi[0] = t.lo[0];
  QVERIFY(!t.isValid());
  QCOMPARE(t.cellCount(), 0LL);
  QCOMPARE(cellOf(t, 0.5, 12.0), -1);

  t = grid2d(0, 5);
  QVERIFY(!t.isValid());
  t = grid2d(4, 5);
  t.bins.pop_back();
  QVERIFY(!t.isValid());
}

void MomeTest::gridCoordsRoundTrip()
{
  const MomeTessellation t = grid2d(4, 5);
  for (int cell = 0; cell < 20; ++cell) {
    const std::vector<int> c = gridCellCoords(cell, t);
    QCOMPARE(c.size(), size_t(2));
    QCOMPARE(c[0] * 5 + c[1], cell);
  }
  QVERIFY(gridCellCoords(20, t).empty());
  QVERIFY(gridCellCoords(-1, t).empty());
}

void MomeTest::gridSupportsOneAndThreeDimensions()
{
  MomeTessellation one;
  one.lo = { 0.0 };
  one.hi = { 10.0 };
  one.bins = { 5 };
  const double v = 7.0;
  QCOMPARE(assignCell(&v, one), 3);

  MomeTessellation three;
  three.lo = { 0.0, 0.0, 0.0 };
  three.hi = { 1.0, 1.0, 1.0 };
  three.bins = { 2, 3, 4 };
  QCOMPARE(three.cellCount(), 24LL);
  const double p[3] = { 0.75, 0.5, 0.9 };
  QCOMPARE(assignCell(p, three), 1 * 12 + 1 * 4 + 3);
  QCOMPARE(gridCellCoords(19, three), std::vector<int>({ 1, 1, 3 }));
}

void MomeTest::cvtCentroidsAreDeterministicForSeed()
{
  const auto a = buildCvtCentroids(2, 16, 2000, 7);
  const auto b = buildCvtCentroids(2, 16, 2000, 7);
  const auto c = buildCvtCentroids(2, 16, 2000, 8);
  QCOMPARE(a.size(), size_t(16));
  QVERIFY(a == b);
  QVERIFY(a != c);

  QVERIFY(buildCvtCentroids(0, 16, 2000, 7).empty());
  QVERIFY(buildCvtCentroids(2, 0, 2000, 7).empty());
  QCOMPARE(buildCvtCentroids(2, 50, 10, 7).size(), size_t(10)); // cells capped by samples
}

void MomeTest::cvtCentroidsSpreadOverUnitBox()
{
  const auto cents = buildCvtCentroids(2, 4, 20000, 1);
  // Four centroids of a uniform square settle near the quadrant centres:
  //   every quadrant holds exactly one.
  std::set<int> quadrants;
  for (const auto& c : cents) {
    QVERIFY(c[0] > 0.0 && c[0] < 1.0 && c[1] > 0.0 && c[1] < 1.0);
    quadrants.insert((c[0] > 0.5 ? 2 : 0) + (c[1] > 0.5 ? 1 : 0));
  }
  QCOMPARE(quadrants.size(), size_t(4));
}

void MomeTest::cvtAssignsNearestCentroid()
{
  MomeTessellation t;
  t.kind = MomeTessellation::Cvt;
  t.lo = { 0.0, 10.0 };
  t.hi = { 1.0, 20.0 };
  t.centroids = { { 0.25, 0.25 }, { 0.75, 0.25 }, { 0.5, 0.9 } };
  QVERIFY(t.isValid());
  QCOMPARE(t.cellCount(), 3LL);
  QCOMPARE(cellOf(t, 0.2, 12.0), 0);  // normalized (0.2, 0.2)
  QCOMPARE(cellOf(t, 0.8, 12.0), 1);
  QCOMPARE(cellOf(t, 0.5, 19.0), 2);
  QCOMPARE(cellOf(t, 5.0, -50.0), 1); // clamped to (1, 0)
  QVERIFY(gridCellCoords(0, t).empty());

  t.centroids.clear();
  QVERIFY(!t.isValid());
  QCOMPARE(cellOf(t, 0.2, 12.0), -1);
}

void MomeTest::truncationKeepsBoundaryAndDropsMostCrowded()
{
  // A 2-objective front; points 2 and 3 are nearly on top of each other.
  const std::vector<std::vector<double>> obj = {
    { 0.0, 10.0 }, { 2.0, 8.0 }, { 4.0, 6.0 }, { 4.1, 5.9 }, { 7.0, 3.0 }, { 10.0, 0.0 }
  };
  std::mt19937 rng(0);
  const std::vector<int> kept = truncateFrontByCrowding(obj, { 0, 1, 2, 3, 4, 5 }, 5, rng);
  QCOMPARE(kept.size(), size_t(5));
  QVERIFY(std::count(kept.begin(), kept.end(), 0) == 1); // boundary
  QVERIFY(std::count(kept.begin(), kept.end(), 5) == 1); // boundary
  // exactly one of the crowded pair is gone
  QCOMPARE(std::count(kept.begin(), kept.end(), 2) + std::count(kept.begin(), kept.end(), 3), 1L);

  // Down to 2: only the boundary points remain.
  const std::vector<int> two = truncateFrontByCrowding(obj, { 0, 1, 2, 3, 4, 5 }, 2, rng);
  QCOMPARE(two, std::vector<int>({ 0, 5 }));
}

void MomeTest::truncationLeavesSmallFrontsAlone()
{
  const std::vector<std::vector<double>> obj = { { 0.0, 1.0 }, { 1.0, 0.0 }, { 0.5, 0.5 } };
  std::mt19937 rng(0);
  QCOMPARE(truncateFrontByCrowding(obj, { 2, 0, 1 }, 3, rng), std::vector<int>({ 2, 0, 1 }));
  QCOMPARE(truncateFrontByCrowding(obj, { 2, 0, 1 }, 10, rng), std::vector<int>({ 2, 0, 1 }));
  // maxSize below 1 is treated as 1
  QCOMPARE(truncateFrontByCrowding(obj, { 2, 0, 1 }, 0, rng).size(), size_t(1));
}

void MomeTest::weightsReplaceInfiniteWithMaxFinite()
{
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<double> w = crowdingSelectionWeights({ inf, 0.5, 1.5, inf });
  QCOMPARE(w, std::vector<double>({ 1.5, 0.5, 1.5, 1.5 }));
}

void MomeTest::weightsUniformForSmallOrDegenerateFronts()
{
  const double inf = std::numeric_limits<double>::infinity();
  QCOMPARE(crowdingSelectionWeights({ inf }), std::vector<double>({ 1.0 }));
  QCOMPARE(crowdingSelectionWeights({ inf, inf }), std::vector<double>({ 1.0, 1.0 }));
  QCOMPARE(crowdingSelectionWeights({ inf, inf, inf }), std::vector<double>({ 1.0, 1.0, 1.0 }));
  QCOMPARE(crowdingSelectionWeights({ 0.0, 0.0, 0.0 }), std::vector<double>({ 1.0, 1.0, 1.0 }));
  QCOMPARE(sampleByCrowding({}, 0.3), -1);
  QCOMPARE(sampleByCrowding({ inf, inf }, 0.0), 0);
  QCOMPARE(sampleByCrowding({ inf, inf }, 0.75), 1);
}

void MomeTest::samplingFrequenciesFollowCrowding()
{
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<double> dists = { inf, 0.5, 1.5, 0.0, inf }; // weights 1.5 .5 1.5 0 1.5
  const std::vector<double> expected = { 0.3, 0.1, 0.3, 0.0, 0.3 };

  std::mt19937 rng(42);
  std::vector<int> hits(dists.size(), 0);
  const int n = 200000;
  for (int i = 0; i < n; ++i) {
    const int k = sampleByCrowding(dists, static_cast<double>(rng()) / 4294967296.0);
    QVERIFY(k >= 0 && k < static_cast<int>(dists.size()));
    ++hits[static_cast<size_t>(k)];
  }
  for (size_t i = 0; i < dists.size(); ++i)
    QVERIFY2(std::fabs(static_cast<double>(hits[i]) / n - expected[i]) < 0.01,
             qPrintable(QString("index %1: %2").arg(i).arg(static_cast<double>(hits[i]) / n)));

  // Edge values of r
  QCOMPARE(sampleByCrowding(dists, 0.0), 0);
  QCOMPARE(sampleByCrowding(dists, 1.0), 4);
}

void MomeTest::archiveKeepsFirstFrontPerCell()
{
  const MomeTessellation t = grid2d(2, 1); // two cells split at x = 0.5
  const std::vector<std::vector<double>> desc = {
    { 0.1, 15.0 }, { 0.2, 15.0 }, { 0.3, 15.0 }, // cell 0
    { 0.9, 15.0 }, { 0.8, 15.0 }                 // cell 1
  };
  const std::vector<std::vector<double>> obj = {
    { 1.0, 5.0 }, { 5.0, 1.0 }, { 6.0, 6.0 }, // row 2 is dominated inside cell 0
    { 9.0, 9.0 }, { 8.0, 8.0 }                // row 3 is dominated by row 4
  };
  const MomeArchive a = buildMomeArchive(obj, desc, t, 10);

  QCOMPARE(a.cellOf, std::vector<int>({ 0, 0, 0, 1, 1 }));
  QCOMPARE(a.cells.size(), size_t(2));
  QCOMPARE(a.cells[0].cell, 0);
  QCOMPARE(a.cells[0].members, std::vector<int>({ 0, 1 }));
  QCOMPARE(a.cells[0].crowding.size(), size_t(2));
  QCOMPARE(a.cells[1].cell, 1);
  // Row 4 would be dominated globally by nothing in its own cell: it stays,
  //   even though rows 0 and 1 are "better" elsewhere in the grid.
  QCOMPARE(a.cells[1].members, std::vector<int>({ 4 }));
}

void MomeTest::archiveRespectsMaxFrontSize()
{
  const MomeTessellation t = grid2d(1, 1);
  std::vector<std::vector<double>> desc, obj;
  for (int i = 0; i <= 20; ++i) {
    desc.push_back({ 0.5, 15.0 });
    obj.push_back({ static_cast<double>(i), static_cast<double>(20 - i) });
  }
  const MomeArchive a = buildMomeArchive(obj, desc, t, 6);
  QCOMPARE(a.cells.size(), size_t(1));
  const MomeCell& c = a.cells[0];
  QCOMPARE(c.members.size(), size_t(6));
  QCOMPARE(c.crowding.size(), size_t(6));
  QCOMPARE(c.members.front(), 0);  // both extremes survive
  QCOMPARE(c.members.back(), 20);
  for (size_t i = 0; i < c.members.size(); ++i) {
    for (size_t j = 0; j < c.members.size(); ++j) {
      if (i != j)
        QVERIFY(!doesDominate(obj[static_cast<size_t>(c.members[i])].data(),
                              obj[static_cast<size_t>(c.members[j])].data(), 2));
    }
  }
  // Same input and seed: same archive.
  const MomeArchive b = buildMomeArchive(obj, desc, t, 6);
  QCOMPARE(b.cells[0].members, c.members);
}

void MomeTest::archiveSkipsUnassignableRows()
{
  const MomeTessellation t = grid2d(2, 2);
  const std::vector<std::vector<double>> desc = {
    { 0.1, 11.0 }, { std::numeric_limits<double>::quiet_NaN(), 11.0 }, { 0.1 }
  };
  const std::vector<std::vector<double>> obj = { { 1.0 }, { 0.0 }, { 0.0 } };
  const MomeArchive a = buildMomeArchive(obj, desc, t, 10);
  QCOMPARE(a.cellOf, std::vector<int>({ 0, -1, -1 }));
  QCOMPARE(a.cells.size(), size_t(1));
  QCOMPARE(a.cells[0].members, std::vector<int>({ 0 }));

  MomeTessellation bad = t;
  bad.bins.clear();
  QVERIFY(buildMomeArchive(obj, desc, bad, 10).cells.empty());
}

QTEST_GUILESS_MAIN(MomeTest)

#include "mometest.moc"
