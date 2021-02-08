//-----------------------------------------------------------------------------
// gmslib - Gaussian Mixture Surface Library
// Copyright (c) Reinhold Preiner 2014-2020 
// 
// Usage is subject to the terms of the WFP (modified BSD-3-Clause) license.
// See the accompanied LICENSE file or
// https://github.com/rpreiner/gmslib/blob/main/LICENSE
//-----------------------------------------------------------------------------

#pragma once

#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>

#include "gaussian.hpp"
#include "pointindex.hpp"
#include "sphereindex.hpp"
#include "pointset.hpp"
#include "random.hpp"
#include "timer.hpp"
#include "memoryinfo.hpp"
#include <omp.h>

using namespace std;



#define GMS_DEBUG	0


namespace gms
{
	/// Mixture of Gaussians
	class Mixture : public vector<Gaussian>
	{
	public:
		// models the distribution of the Gaussian's normals using a spherical Gaussian
		// direction: mean normal. length: variance of spherical Gaussian + 1 
		// (the +1 offset allows to encode the variance in the normal length, while representing zero variance as valid vector)
		vector<vec3> nvars;

		bool hasNormals() const { return !nvars.empty(); }

		using ptr = shared_ptr<Mixture>;

	public:
		struct Params
		{
			bool	verbose = false;				// verbose console output
			bool	memoryProfiling = false;		// activate memory usage profiling
			int		initNeighborhoodType = 1;		// 0: initialize Gaussian using all samples within maxInitNeighborDist, 1: use only kNNCount nearest neighbors.
			uint	kNNCount = 8;					// number nearest neighbors per point used for initial Gaussian computation. The neighbor set is clamped by maxInitNeighborRadius.
			float	maxInitNeighborDist = 1.0f;		// global initialization kernel radius
			float	initIsotropicStdev = 1.0f;		
			bool	initIsotropic = false;			// isotropic initial Gaussians, with stddev initIsotropicStdev
			//bool	useGlobalInitRadius = true;		// use global initialization radius instead of NNDist sampling
			//uint	nNNDistSamples = 10;			// number of sample points for computing the grid cell size based on nearest neighbor distances
			bool	useWeightedPotentials = true;	// if true, performs WLOP-like balancing of the initial Gaussian potentials
			bool	initMeansInPoints = true;		// positions the initial Gaussians in the point positions instead of the local means
			uint	nLevels = 4;					// number of levels to use when clustering
			float	hemReductionFactor = 3.0f;		// factor by which to reduce the mixture each level
			float	alpha = 2.2f;					// multiple of cluster maximum std deviation to use for query radius
			bool	computeNVar = true;
			bool	blockProcessing = false;
			uint	blockSize = 1000000;
			uint	numThreads;
			Params() {}
		};
				
	public:
		// Default Constructor
		Mixture() {}


		// Constructor. Takes a set of 3D points and computes a mixture using HEM
		Mixture(const PointSet* points, const Params& params = Params()) {
			create(points, params);
		}

		Mixture(const Mixture& M) {
			*this = M;
		}

		Mixture(uint nGaussians) {
			resize(nGaussians);
		}


		const Mixture& operator= (const Mixture& M)
		{
			if (&M == this)
				return *this;

			// copy components
			clear();
			for (const Gaussian& g : M)
				push_back(g);

			// copy attributes
			if (M.hasNormals())	nvars = M.nvars;
						
			Memory::instance()->record(__LINE__);

			return *this;
		}


		template<typename booltype>
		void removeElements(vector<booltype>& invalidGaussians)
		{
			unordered_map<uint, uint> gidMap, pidMap;

			// copying necessary data into temporary mixture and swapping back
			{
				Mixture newMixture;

				// copy valid Gaussians
				for (uint s = 0; s < size(); s++)
				{
					if (!invalidGaussians[s])
					{
						uint gid = newMixture.size();
						newMixture.push_back(at(s));
						if (hasNormals())	newMixture.nvars.push_back(nvars[s]);
						gidMap[s] = gid;
					}
				}
				
				Memory::instance()->record(__LINE__);

				// at this point we can swap back the compactified mixture info - old one is deallocated
				this->swap(newMixture);
				nvars.swap(newMixture.nvars);
			}
			Memory::instance()->record(__LINE__);
		}


		void create(const PointSet* points, const Params& params = Params())
		{
			random::reset();
			
			// set number of threads used
			omp_set_dynamic(0);
			omp_set_num_threads(params.numThreads);
			#pragma omp parallel
			#pragma omp master
			{ cout << "using " << omp_get_num_threads() << " threads" << endl; }

			if (params.memoryProfiling)
				Memory::instance()->setActive(true);

			// Blocked Processing: build a kd-tree on the points and process leaf by leaf
			if (params.blockProcessing && params.blockSize > 0)
			{
				vector<const vector<uint>*> leaves;
				SphereIndex kdTree(*points, 0, params.blockSize);
				kdTree.getLeafIndices(leaves);

				cout << "Created kd-tree with " << leaves.size() << " blocks." << endl << "Building mixtures ";

				// build mixture leaf by leaf
				clear();
				nvars.clear();
				
				for (uint i = 0; i < leaves.size(); i++)
				{
					cout << "---- Mixture #" << i << " ----" << endl;
					const vector<uint>* leaf = leaves[i];
										
					PointSet leafPoints;
					for (const string& name : points->getFloatAttributeNames())	leafPoints.setAttrib(name, make_shared<vector<float>>());
					for (const string& name : points->getVec3AttributeNames())	leafPoints.setAttrib(name, make_shared<vector<gms::vec3>>());
					for (const string& name : points->getVec4AttributeNames())	leafPoints.setAttrib(name, make_shared<vector<gms::vec4>>());

					for (uint j : *leaf)
					{
						leafPoints.push_back(points->at(j));

						for (auto& attrib : points->getFloatAttribs())	leafPoints.attribFloat(attrib.first)->push_back(attrib.second->at(j));
						for (auto& attrib : points->getVec3Attribs())	leafPoints.attribVec3(attrib.first)->push_back(attrib.second->at(j));
						for (auto& attrib : points->getVec4Attribs())	leafPoints.attribVec4(attrib.first)->push_back(attrib.second->at(j));
					}

					// create mixture
					Mixture leafMixture;
					leafMixture.createBlock(&leafPoints, params);

					// append mixture to this mixture
					insert(end(), leafMixture.begin(), leafMixture.end());
					if (leafMixture.hasNormals())	nvars.insert(nvars.end(), leafMixture.nvars.begin(), leafMixture.nvars.end());
				}
			}
			// create mixture from whole point set
			else
			{
				createBlock(points, params);
			}
		}


		// reduces the mixture by nLevels levels. If nLevels is zero, the mixture is reduced until no further reduction is possible.
		void reduce(const Params& params)
		{
			//cout << "Reducing mixture with size " << size() << endl;
			Timer timer; 
			timer.start();
			// complete reduction until convergence
			if (params.nLevels == 0)
			{
				uint lastSize = 0;
				uint l = 1;
				while (size() != lastSize)
				{
					if (params.verbose)	cout << "level " << l;
					++l;
					lastSize = size();
					reduceLevel(params);
					if (params.verbose) cout << ":\tsize " << size() << endl;
				}
				cout << (l - 1) << " levels\n";
			}
			// reduction by n levels
			else
			{
				for (uint l = 1; l <= params.nLevels; ++l)
				{
					if (params.verbose)	cout << "level " << l;
					reduceLevel(params);
					if (params.verbose)	cout << ":\tsize " << size() << endl;
				}
			}
			cout << "  Reduced in " << timer.stop() << " ms" << endl;
		}



		struct SelectChildrenProc : PointIndex::PointProcessor
		{
			// params
			const Mixture& mixture;
			const vector<char> isParent;
			const vector<uint>& packedParentPos;
			const vector<float>& queryRadii;
			uint parentIndex;
			float half_alpha2;
			// output
			vector<vector<uint>>& outChildIndices;


			SelectChildrenProc(
				const Mixture& mixture, 
				const vector<char>& isParent, 
				const vector<uint>& packedParentPos, 
				const vector<float>& queryRadii, 
				float half_alpha2, 
				vector<vector<uint>>& outChildIndices
			) : 
			mixture(mixture), 
			isParent(isParent), 
			packedParentPos(packedParentPos), 
			queryRadii(queryRadii), 
			half_alpha2(half_alpha2), 
			outChildIndices(outChildIndices)
			{}
			
			virtual void operator() (uint index, const vec3& pos, const vector<uint>& neighbors)
			{
				const Gaussian& parent = mixture[index];
				uint packedPos = packedParentPos[index];
				
				// only process parents
				if (!isParent[index])
					return;
				
				// select eligible children from the conservative resultSet
				const float sqRadius = queryRadii[packedPos] * queryRadii[packedPos];
				for (uint i : neighbors)
				{
					const Gaussian& child = mixture[i];
					if (sqdist(parent.mu, child.mu) < sqRadius)
					{

						// consider not taking this child i into account only if its not the parent s itself!
						if (i != index)
						{
							// don't merge other parents or childs past the KLD threshold
							if (isParent[i] || KLD(child, parent) > half_alpha2)
								continue;
						}
				
						outChildIndices[packedPos].push_back(i);
					}
				}
				
				// if the result set is empty due to numerical instabilities at small query radii, ensure that the parent itself is part of the child set!
				// this is necessary to avoid zero sum likelihoods and consequently nan responsibilities!
				if (outChildIndices[packedPos].empty())
				{
					outChildIndices[packedPos].push_back(index);
					return;
				}
			}
		};



		// reduce mixture by one level given the regularization constraint parameter alpha. 
		// The hemReductionFactor is recommended to be 3.0.
		void reduceLevel(const Params& params)
		{
			//Timer timer;
			const float parentProbability = 1.0f / params.hemReductionFactor;
			uint nGaussians = size();

			Memory::instance()->record(__LINE__);

			// 1. iterate over components, prepare index point centers and compute the parent's individual and maximum query radius
			vector<char> isParent(nGaussians, 0);
			vector<vec3> centers(nGaussians);
			vector<uint> packedParentPos(nGaussians, UINT_MAX);
			vector<uint> parentIndices;
			vector<float> queryRadii;
			float maxQueryRadius = 0;

			//timer.start();
			for (uint i = 0; i < nGaussians; ++i)
			{
				Gaussian& g = at(i);

				// prepare centers to be indexed
				centers[i] = g.mu;

				// select parents and save parent flag
				if (isParent[i] = random::uniform01() < parentProbability ? 1 : 0)
				{
					packedParentPos[i] = parentIndices.size();
					parentIndices.push_back(i);

					// ii. get the conservative query radius for this parent
					float queryRadius = params.alpha * sqrtf(g.cov.eigenvalues().z);
					queryRadii.push_back(queryRadius);

					// ii. determine maximum query radius
					if (queryRadius > maxQueryRadius)
						maxQueryRadius = queryRadius;
				}
			}
			uint nParents = parentIndices.size();
			//cout << "init parents set: " << timer.stop() << " ms" << endl;

			Memory::instance()->record(__LINE__);

			// 2. create point index of component centers for neighbor queries
			//timer.start();
			PointIndex* index = new PointIndex(centers, maxQueryRadius);
			//cout << "create index: " << timer.stop() << " ms" << endl;

			Memory::instance()->record(__LINE__);

			// 3. select child set for each parent
			//timer.start();
			vector<vector<uint>> childIndices(nParents);
			vector<vec3i> cellCoords = index->getCellCoords();

			const float half_alpha2 = params.alpha * params.alpha * 0.5f;
			SelectChildrenProc proc(*this, isParent, packedParentPos, queryRadii, half_alpha2, childIndices);
			#pragma omp parallel for
			for (int i = 0; i < (int)cellCoords.size(); i++)
				index->processCell(cellCoords[i], proc);


			//cout << "select child set: " << timer.stop() << " ms" << endl;
			Memory::instance()->record(__LINE__);

			delete index;
			packedParentPos = vector<uint>();
			queryRadii = vector<float>();
			centers = vector<vec3>();
			cellCoords = vector<vec3i>();
			Memory::instance()->record(__LINE__);


			// 4. compute the wL_is and the wL sums
			//timer.start();
			vector<vector<float>> wL_cache(nParents);
			vector<float> sumLw(nGaussians, 0);
			for (int s_ = 0; s_ < (int)nParents; ++s_)
			{
				uint s = parentIndices[s_];
				const Gaussian& parent = at(s);
				
				// iterate over children 
				const vector<uint>& I = childIndices[s_];
				wL_cache[s_].resize(I.size(), 0.0f);
				for (uint i_ = 0; i_ < I.size(); ++i_)
				{
					uint i = I[i_];
					const Gaussian& child = at(i);

					float likelihood = hemLikelihood(parent, child);

					const float maxL = 1e8f;
					const float minL = FLT_MIN;
					float wL_si = parent.weight * clamp(likelihood, minL, maxL);
					
					// save likelihood contribution
					wL_cache[s_][i_] = wL_si;
					
					sumLw[i] += wL_si;
				}
			}
			//cout << "create wlCache: " << timer.stop() << " ms" << endl;

			Memory::instance()->record(__LINE__);


			// 5. compute responsibilities and update
			//timer.start();
			Mixture newMixture(nParents);
			if (hasNormals()) newMixture.nvars.resize(nParents);
			
			#pragma omp parallel for
			for (int s_ = 0; s_ < (int)nParents; ++s_)
			{
				uint s = parentIndices[s_];
				const Gaussian& parent = at(s);
				const vector<uint>& I = childIndices[s_];

				// initialize parent info
				float w_s = 0.0f;
				vec3 summu_i(0, 0, 0);
				smat3 sumcov_i(0, 0, 0, 0, 0, 0);
				vec3 resultant(0, 0, 0);
				float nvar = 0.0f;
				
				// iterate over children and accumulate
				for (uint i_ = 0; i_ < I.size(); ++i_)
				{
					uint i = I[i_];

					if (sumLw[i] == 0.0f)	// can happen
						continue;

					const Gaussian& child = at(i);

					// compute responsibility of parent s for child i
					float r_is = wL_cache[s_][i_] / sumLw[i];
					float w = r_is * child.weight;

					// accumulate
					w_s += w;
					summu_i += w * child.mu;
					sumcov_i += w * (child.cov + smat3::outer(child.mu - parent.mu));	// accumulates generic cov relative to parent mu, numerically more stable than origin, due to smaller distances

					if (hasNormals())
					{
						// normal cluster update
						float c_nvar = length(nvars[i]) - 1;		// decode by substracting the +1 offset
						vec3 c_normal = normalize(nvars[i]);

						// flip child normal to be next to the parent normal (parent.nvar is unnormalized, but thats ok, we're just using its direction)
						if (dot(c_normal, nvars[s]) < 0.0f)
							c_normal = -c_normal;

						resultant += w * c_normal;
						nvar += w * c_nvar;
					}
				}

				// normalize and condition new cov matrix
				float inv_w = 1.0f / w_s;		// w_s > 0 is certain
				vec3 mu_s = inv_w * summu_i;
				smat3 cov_s = inv_w * sumcov_i - smat3::outer(mu_s - parent.mu);
				cov_s = conditionCov(cov_s);

				// Add new component to output list
				Gaussian newComponent;
				newComponent.mu = mu_s;
				newComponent.cov = cov_s;
				newComponent.weight = w_s;
				newMixture[s_] = newComponent;

				// mixture of normals
				if (hasNormals())
				{
					float variance1 = nvar * inv_w;			// normalized sum of the variances of the child clusters
					float R = length(resultant);			// resultant length
					float Rmean = R * inv_w;				// mean resultant length
					float variance2 = -2.0f * log(Rmean);	// variance of the child clusters' mean normal vectors with respect to the common mean (outer product thingy)
					vec3  newMeanNormal = resultant / R;	// normalized mean normal vector of new cluster
					newMixture.nvars[s_] = newMeanNormal * (variance1 + variance2 + 1);	// encode again with +1 offset
				}
			}
			//cout << "update: " << timer.stop() << " ms" << endl;

			Memory::instance()->record(__LINE__);

			childIndices = vector<vector<uint>>();
			wL_cache = vector<vector<float>>();
			Memory::instance()->record(__LINE__);



			// 6. add orphans, components not addressed by any parent (zero sumLw)
			//timer.start();
			vector<int> isOrphan(nGaussians, 0);
			vector<uint> indices(nGaussians);
			vector<uint> packedIndices(nGaussians);
			#pragma omp parallel for
			for (int i = 0; i < (int)nGaussians; ++i)
			{
				indices[i] = i;
				if (sumLw[i] == 0.0f)
					isOrphan[i] = 1;
			}

			uint packedSize = parallel::pack(indices, isOrphan, packedIndices);
			uint oldSize = newMixture.size();
			newMixture.resize(oldSize + packedSize);
			if (hasNormals())	newMixture.nvars.resize(oldSize + packedSize);
			#pragma omp parallel for
			for (int i = 0; i < (int)packedSize; ++i)
			{
				newMixture[oldSize + i] = at(packedIndices[i]);
				if (hasNormals())	newMixture.nvars[oldSize + i] = nvars[packedIndices[i]];
			}
			//cout << "orphans: " << timer.stop() << " ms" << endl;

			Memory::instance()->record(__LINE__);
						
			// 7. exchange current set of components with new reduced one
			this->swap(newMixture);
			if (hasNormals()) nvars.swap(newMixture.nvars);

			Memory::instance()->record(__LINE__);

#if GMS_DEBUG
			for (uint i = 0; i < newMixture.size(); ++i)
			{
				const vec3& mu = newMixture[i].mu;
				const smat3& cov = newMixture[i].cov;
				if (isnan(mu) || det(cov) <= 0 || isnan(det(cov)))
					cerr << "[reduceLevel] Error @ new component " << i << ": mu = " << mu << ", cov = " << cov << ", det = " << det(cov) << endl;
			}
#endif
		}


		static float hemLikelihood(const Gaussian& parent, const Gaussian& child)
		{
			vec3 mu_diff = parent.mu - child.mu;
			smat3 pCovInv = inverse(parent.cov);

			float smd = dot(mu_diff, pCovInv * mu_diff);
			mat3 ipcCov = pCovInv * child.cov;

			// Gaussian exponent and the trace exponent
			float e = -0.5f * (smd + trace(ipcCov));
			// 1/((2pi)^(3/2))   *   sqrt(|Sigma1|)^-1 = sqrt(|Sigma1^-1|)   *   e^exponent   
			float f = 0.063493635934241f * sqrtf(det(pCovInv)) * expf(e);
			// raise to the power of the number of points in the child cluster
			return powf(f, child.weight);
		}
		

	protected:
		// creates a mixture from a whole PointSet block points
		void createBlock(const PointSet* points, const Params& params = Params())
		{
			clear();
			nvars.clear();

			// 1. initialize mixture M
			cout << "initializing ";
			Timer timer;
			timer.start();
			initMixture(*points, params);
			cout << " ... " << timer.stop() << " ms" << endl;

			Memory::instance()->record(__LINE__);

			// 2. hierarchical clustering
			if (params.nLevels > 0)
				reduce(params);

			Memory::instance()->record(__LINE__);

			if (Memory::instance()->isActive())
				cerr << "  Peak Relative Mem Usage: " << Memory::instance()->peakRelativeUsage() << " MB" << endl;
		}

		
		void initMixture(const PointSet& points, const Params& params)
		{
			smat3 initIsotropicCov = (params.initIsotropicStdev * params.initIsotropicStdev) * smat3::identity();

			// in case the initial Gaussians should be isotropic with constant variance and initialized in the point's pos, no need to query neighbors
			if (params.initIsotropic && params.initMeansInPoints)
			{
				resize(points.size());
				#pragma omp	parallel for
				for (int i = 0; i < (int)points.size(); ++i)
					at(i) = Gaussian(points[i], initIsotropicCov);
			}
			else
			{
				bool useGlobalRadius = params.initNeighborhoodType == 0;
				PointIndex index(points, params.maxInitNeighborDist);

				// Compute initial mixture
				resize(points.size());

#pragma omp	parallel for
				for (int i = 0; i < (int)points.size(); ++i)
				{
					const vec3& q = points[i];
					Gaussian& g = at(i);


					// query neighbors using knn search, constrained by maxInitNeighborDist
					vector<uint> nnIndices;
					if (useGlobalRadius)
						index.radiusSearch(q, params.maxInitNeighborDist, nnIndices);
					else
						index.annSearch(q, max(params.kNNCount, 1u), nnIndices);

					// in case we didn't find any neighbor
					if (nnIndices.size() < 2)
					{
						g = Gaussian(q, smat3::identity() * params.maxInitNeighborDist * 0.0001f);
						continue;
					}

					// compute average squared distance to all neighbors
					float inv_w = 1.0f / nnIndices.size();		// size > 0 ensured
					float avgdist = 0;
					for (uint j : nnIndices)
						avgdist += sqdist(q, points[j]);
					avgdist *= inv_w;

					// compute cov
					const float minus_16_over_r2 = -16.0f / (params.maxInitNeighborDist * params.maxInitNeighborDist);
					float eps = avgdist * avgdist * 0.0001f;
					smat3 sumCov = smat3(eps, 0, 0, eps, 0, eps);	// initial bias of the cluster for stability
					vec3 sumMean(0, 0, 0);
					float density = 0.000001f;
					for (uint j : nnIndices)
					{
						const vec3& p = points[j];
						vec3 diff = p - q;
						sumCov += smat3::outer(diff);				// centralized in parent pos for numerical stability
						sumMean += p;
						if (useGlobalRadius && params.useWeightedPotentials)
							density += expf(minus_16_over_r2 * dot(diff, diff)); 	// accumulate the WLOP density sum
					}

					// setup component
					vec3 mean = sumMean * inv_w;
					g.mu = params.initMeansInPoints ? q : mean;
					g.weight = useGlobalRadius && params.useWeightedPotentials ? 1.0f / density : 1.0f;	

					// in case the initial Gaussians should be isotropic with constant variance
					g.cov = initIsotropicCov;

					if (!params.initIsotropic)
					{
						g.cov += sumCov * inv_w - smat3::outer(mean - q);	// consider parent pos centralization
						g.cov = conditionCov(g.cov);
					}
				}
			}

			// initialize distributions of normals
			if (params.computeNVar)
			{
				cout << "initializing normals" << endl;
				const float initialVar = 0.001f;	// for the initial normal variance, 0.0f should theoretically also work

				nvars.resize(points.size());
				#pragma omp parallel for
				for (int i = 0; i < (int)points.size(); i++)
				{
					// Compute the initial normal and set initial normal variance of this point cluster
					// the normal variance is encoded in the length of the normal vector
					vec3 evectors[3];
					at(i).cov.eigenvectors(evectors);
					nvars[i] = evectors[0] * (initialVar + 1);
				}
			}
				
		#if GMS_DEBUG
			for (uint i = 0; i < size(); ++i)
			{
				const vec3& mu = at(i).mu;
				const smat3& cov = at(i).cov;
				if (isnan(mu) || det(cov) <= 0 || isnan(det(cov)))
					cerr << "[initMixture] Error @ initial component " << i << ": mu = " << mu << ", cov = " << cov << ", det = " << det(cov) << endl;
			}
		#endif
		}



#pragma region operations
	public:
		// returns the mixture pdf at the given point x using the Gaussians of the given indices
		float pdf(vector<uint>& indices, const vec3& x) const
		{
			float pdf = 0.0f;
			for (uint s : indices)
			{
				const Gaussian& g = at(s);
				pdf += g.weight * g.pdf(x);
			}
			return pdf;
		}

		// returns the gradient of the mixture pdf at the given point x using the Gaussians of the given indices
		vec3 gradient(vector<uint>& indices, const vec3& x) const
		{
			vec3 grad(0, 0, 0);
			for (uint s : indices)
			{
				const Gaussian& g = at(s);
				grad += g.weight * g.gradient(x);
			}
			return grad;
		}
#pragma endregion

	};	/// end class Mixture


}	/// end namespace gms

