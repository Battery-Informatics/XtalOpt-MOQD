/**********************************************************************
  mome - MOME-X (crowding-based Multi-Objective MAP-Elites) archive math:
         descriptor-space tessellation (grid or CVT), per-cell size-limited
         Pareto fronts, and crowding-biased sampling.

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

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace Search {

namespace {

// Uniform double in [0, 1), identical on every platform for a given engine state.
double uniform01(std::mt19937& rng)
{
  return static_cast<double>(rng()) / 4294967296.0;
}

// Map a raw descriptor value to [0, 1] within its bounds (clamped).
double normalizedClamped(double v, double lo, double hi)
{
  const double t = (v - lo) / (hi - lo);
  return t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
}

int nearestCentroid(const double* p, const std::vector<std::vector<double>>& centroids, int dims)
{
  int best = -1;
  double bestDist = std::numeric_limits<double>::infinity();
  for (size_t c = 0; c < centroids.size(); ++c) {
    double d = 0.0;
    for (int k = 0; k < dims; ++k) {
      const double diff = p[k] - centroids[c][static_cast<size_t>(k)];
      d += diff * diff;
    }
    if (d < bestDist) {
      bestDist = d;
      best = static_cast<int>(c);
    }
  }
  return best;
}

} // namespace

long long MomeTessellation::cellCount() const
{
  if (!isValid())
    return 0;
  if (kind == Cvt)
    return static_cast<long long>(centroids.size());
  long long n = 1;
  for (int b : bins)
    n *= b;
  return n;
}

bool MomeTessellation::isValid() const
{
  const size_t d = lo.size();
  if (d == 0 || hi.size() != d)
    return false;
  for (size_t i = 0; i < d; ++i) {
    if (!std::isfinite(lo[i]) || !std::isfinite(hi[i]) || !(lo[i] < hi[i]))
      return false;
  }
  if (kind == Grid) {
    if (bins.size() != d)
      return false;
    long long n = 1;
    for (int b : bins) {
      if (b < 1)
        return false;
      n *= b;
      if (n > std::numeric_limits<int>::max())
        return false;
    }
    return true;
  }
  if (centroids.empty())
    return false;
  for (const auto& c : centroids) {
    if (c.size() != d)
      return false;
  }
  return true;
}

std::vector<std::vector<double>> buildCvtCentroids(int dims, int cells, int samples,
                                                   unsigned int seed, int iterations)
{
  if (dims < 1 || cells < 1 || samples < 1)
    return {};
  if (cells > samples)
    cells = samples;

  std::mt19937 rng(seed);
  const size_t d = static_cast<size_t>(dims);
  std::vector<double> pts(static_cast<size_t>(samples) * d);
  for (double& v : pts)
    v = uniform01(rng);

  // The samples are i.i.d. uniform, so the first ones are a random initial guess.
  std::vector<std::vector<double>> centroids(static_cast<size_t>(cells), std::vector<double>(d));
  for (size_t c = 0; c < centroids.size(); ++c)
    std::copy(pts.begin() + static_cast<long>(c * d), pts.begin() + static_cast<long>((c + 1) * d),
              centroids[c].begin());

  std::vector<double> sums(centroids.size() * d);
  std::vector<int> counts(centroids.size());
  for (int it = 0; it < iterations; ++it) {
    std::fill(sums.begin(), sums.end(), 0.0);
    std::fill(counts.begin(), counts.end(), 0);
    for (int s = 0; s < samples; ++s) {
      const double* p = &pts[static_cast<size_t>(s) * d];
      const size_t c = static_cast<size_t>(nearestCentroid(p, centroids, dims));
      for (size_t k = 0; k < d; ++k)
        sums[c * d + k] += p[k];
      ++counts[c];
    }
    for (size_t c = 0; c < centroids.size(); ++c) {
      if (counts[c] == 0)
        continue; // keep a centroid that attracted no sample where it is
      for (size_t k = 0; k < d; ++k)
        centroids[c][k] = sums[c * d + k] / counts[c];
    }
  }
  return centroids;
}

int assignCell(const double* desc, const MomeTessellation& tess)
{
  if (!desc || !tess.isValid())
    return -1;
  const int d = tess.dims();
  for (int i = 0; i < d; ++i) {
    if (!std::isfinite(desc[i]))
      return -1;
  }

  if (tess.kind == MomeTessellation::Grid) {
    long long cell = 0;
    for (int i = 0; i < d; ++i) {
      const size_t k = static_cast<size_t>(i);
      const int nb = tess.bins[k];
      int b = static_cast<int>(std::floor(normalizedClamped(desc[i], tess.lo[k], tess.hi[k]) * nb));
      b = std::min(std::max(b, 0), nb - 1);
      cell = cell * nb + b;
    }
    return static_cast<int>(cell);
  }

  std::vector<double> p(static_cast<size_t>(d));
  for (int i = 0; i < d; ++i) {
    const size_t k = static_cast<size_t>(i);
    p[k] = normalizedClamped(desc[i], tess.lo[k], tess.hi[k]);
  }
  return nearestCentroid(p.data(), tess.centroids, d);
}

std::vector<int> gridCellCoords(int cell, const MomeTessellation& tess)
{
  if (tess.kind != MomeTessellation::Grid || !tess.isValid() || cell < 0 ||
      cell >= tess.cellCount())
    return {};
  std::vector<int> coords(tess.bins.size());
  for (size_t i = tess.bins.size(); i-- > 0;) {
    coords[i] = cell % tess.bins[i];
    cell /= tess.bins[i];
  }
  return coords;
}

std::vector<int> truncateFrontByCrowding(const std::vector<std::vector<double>>& objData,
                                         std::vector<int> members, int maxSize,
                                         std::mt19937& rng)
{
  if (maxSize < 1)
    maxSize = 1;

  while (static_cast<int>(members.size()) > maxSize) {
    std::vector<std::vector<double>> front;
    front.reserve(members.size());
    for (int m : members)
      front.push_back(objData[static_cast<size_t>(m)]);
    const std::vector<double> dist = frontCrowdingDistance(front);

    const double smallest = *std::min_element(dist.begin(), dist.end());
    std::vector<size_t> ties;
    for (size_t i = 0; i < dist.size(); ++i) {
      if (dist[i] == smallest)
        ties.push_back(i);
    }
    const size_t drop = ties[ties.size() == 1 ? 0 : rng() % ties.size()];
    members.erase(members.begin() + static_cast<long>(drop));
  }
  return members;
}

std::vector<double> crowdingSelectionWeights(const std::vector<double>& dists)
{
  std::vector<double> w(dists.size(), 1.0);
  if (dists.size() <= 2)
    return w;

  double maxFinite = 0.0;
  for (double d : dists) {
    if (std::isfinite(d) && d > maxFinite)
      maxFinite = d;
  }
  if (!(maxFinite > 0.0))
    return w;

  for (size_t i = 0; i < dists.size(); ++i) {
    if (std::isfinite(dists[i]))
      w[i] = dists[i] > 0.0 ? dists[i] : 0.0;
    else
      w[i] = maxFinite; // boundary point (or NaN): as good as the most isolated interior point
  }
  return w;
}

int sampleByCrowding(const std::vector<double>& dists, double r)
{
  if (dists.empty())
    return -1;
  const std::vector<double> w = crowdingSelectionWeights(dists);
  double total = 0.0;
  for (double v : w)
    total += v;

  const double target = r * total;
  double acc = 0.0;
  int lastPositive = 0;
  for (size_t i = 0; i < w.size(); ++i) {
    if (w[i] <= 0.0)
      continue;
    acc += w[i];
    lastPositive = static_cast<int>(i);
    if (target < acc)
      return lastPositive;
  }
  return lastPositive; // r == 1 or rounding
}

MomeArchive buildMomeArchive(const std::vector<std::vector<double>>& objData,
                             const std::vector<std::vector<double>>& descData,
                             const MomeTessellation& tess, int maxSize, unsigned int seed)
{
  MomeArchive archive;
  const size_t n = std::min(objData.size(), descData.size());
  archive.cellOf.assign(n, -1);
  if (!tess.isValid())
    return archive;

  const size_t d = static_cast<size_t>(tess.dims());
  std::map<int, std::vector<int>> byCell; // ordered, so the archive is deterministic
  for (size_t i = 0; i < n; ++i) {
    if (descData[i].size() != d || objData[i].empty())
      continue;
    const int cell = assignCell(descData[i].data(), tess);
    archive.cellOf[i] = cell;
    if (cell >= 0)
      byCell[cell].push_back(static_cast<int>(i));
  }

  std::mt19937 rng(seed);
  for (const auto& entry : byCell) {
    const std::vector<int>& candidates = entry.second;

    std::vector<std::vector<double>> points;
    points.reserve(candidates.size());
    for (int c : candidates)
      points.push_back(objData[static_cast<size_t>(c)]);

    const std::vector<std::vector<int>> fronts = nonDominatedSorting(points);
    if (fronts.empty() || fronts[0].empty())
      continue;

    std::vector<int> local = fronts[0];
    std::sort(local.begin(), local.end());
    std::vector<int> members;
    members.reserve(local.size());
    for (int l : local)
      members.push_back(candidates[static_cast<size_t>(l)]);

    MomeCell cell;
    cell.cell = entry.first;
    cell.members = truncateFrontByCrowding(objData, members, maxSize, rng);

    std::vector<std::vector<double>> front;
    front.reserve(cell.members.size());
    for (int m : cell.members)
      front.push_back(objData[static_cast<size_t>(m)]);
    cell.crowding = frontCrowdingDistance(front);

    archive.cells.push_back(std::move(cell));
  }
  return archive;
}

} // namespace Search
