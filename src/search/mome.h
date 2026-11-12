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

#ifndef SEARCH_MOME_H
#define SEARCH_MOME_H

#include <random>
#include <vector>

namespace Search {

/**
 * Tessellation of the descriptor space into cells.
 *
 * Grid: dimension i is split into bins[i] equal bins over [lo[i], hi[i]];
 *   cells are numbered row-major (the last dimension varies fastest).
 * Cvt: cells are the Voronoi regions of @c centroids, which live in the
 *   bounds-normalized space [0,1]^d (see buildCvtCentroids()).
 *
 * Descriptor values outside [lo, hi] are clamped onto the boundary, so
 * every finite descriptor vector belongs to exactly one cell.
 */
struct MomeTessellation
{
  enum Kind { Grid, Cvt };

  Kind kind = Grid;
  std::vector<double> lo;
  std::vector<double> hi;
  std::vector<int> bins;
  std::vector<std::vector<double>> centroids;

  /** @return Number of descriptor dimensions. */
  int dims() const { return static_cast<int>(lo.size()); }

  /** @return Total number of cells (0 if the tessellation is not usable). */
  long long cellCount() const;

  /** @return Whether bounds, bins/centroids are consistent and usable. */
  bool isValid() const;
};

/**
 * Build @p cells CVT centroids in [0,1]^@p dims with Lloyd's k-means over
 * @p samples uniform sample points. Fully determined by @p seed (the same
 * seed gives the same centroids on every platform).
 *
 * @p cells is reduced to @p samples if larger. Returns an empty list for
 * non-positive @p dims, @p cells, or @p samples.
 */
std::vector<std::vector<double>> buildCvtCentroids(int dims, int cells, int samples,
                                                   unsigned int seed, int iterations = 20);

/**
 * @return The cell index of descriptor vector @p desc (of tess.dims()
 * values), or -1 if the tessellation is invalid or any value is not finite.
 */
int assignCell(const double* desc, const MomeTessellation& tess);

/** @return Per-dimension bin indices of grid cell @p cell (empty for CVT/invalid). */
std::vector<int> gridCellCoords(int cell, const MomeTessellation& tess);

/**
 * Reduce a set of mutually non-dominated points to at most @p maxSize.
 *
 * @p members are row indices into @p objData. While there are too many,
 * the member with the smallest crowding distance (see
 * frontCrowdingDistance()) is dropped and distances are recomputed; ties
 * are broken at random with @p rng. Boundary points (infinite distance)
 * are therefore dropped last.
 *
 * @return The surviving members, in their original relative order.
 */
std::vector<int> truncateFrontByCrowding(const std::vector<std::vector<double>>& objData,
                                         std::vector<int> members, int maxSize,
                                         std::mt19937& rng);

/**
 * Turn raw crowding distances of one front into selection weights:
 * infinite (boundary) distances are replaced by the largest finite
 * distance of the front. If the front has at most 2 members, or no
 * distance is finite and positive, all weights are 1 (uniform).
 */
std::vector<double> crowdingSelectionWeights(const std::vector<double>& dists);

/**
 * Pick an index with probability proportional to the weights given by
 * crowdingSelectionWeights(@p dists), using the uniform random number
 * @p r in [0, 1). Returns -1 for an empty list.
 */
int sampleByCrowding(const std::vector<double>& dists, double r);

/** One occupied cell of the archive. */
struct MomeCell
{
  int cell = -1;                // cell index in the tessellation
  std::vector<int> members;     // row indices into the objective/descriptor data
  std::vector<double> crowding; // raw crowding distance of each member in this front
};

struct MomeArchive
{
  std::vector<MomeCell> cells; // occupied cells only, ordered by cell index
  std::vector<int> cellOf;     // cell index of every input row (-1: not assignable)
};

/**
 * Build the MOME archive from scratch.
 *
 * Row i of @p objData holds the objective values of candidate i (lower is
 * better in every objective), row i of @p descData its raw descriptor
 * values. Each candidate is assigned to a cell; every occupied cell keeps
 * the first Pareto front of its candidates, truncated to @p maxSize by
 * crowding (see truncateFrontByCrowding(); @p seed drives its tie-breaks).
 */
MomeArchive buildMomeArchive(const std::vector<std::vector<double>>& objData,
                             const std::vector<std::vector<double>>& descData,
                             const MomeTessellation& tess, int maxSize,
                             unsigned int seed = 0);

} // namespace Search

#endif // SEARCH_MOME_H
