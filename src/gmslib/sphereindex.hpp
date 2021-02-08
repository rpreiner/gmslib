//-----------------------------------------------------------------------------
// gmslib - Gaussian Mixture Surface Library
// Copyright (c) Reinhold Preiner 2014-2020 
// 
// Usage is subject to the terms of the WFP (modified BSD-3-Clause) license.
// See the accompanied LICENSE file or
// https://github.com/rpreiner/gmslib/blob/main/LICENSE
//-----------------------------------------------------------------------------


#pragma once

#include "vec.hpp"
#include <vector>
#include <algorithm>
#include <cassert>
#include "geom.hpp"
#include "parallel.hpp"

using namespace std;


namespace gms
{
	
	// kd-Tree poluplated by bounding spheres
	class SphereIndex
	{
	private:
		const static uint DEFAULT_LEAF_SIZE = 2048;		// empirically determined to be somewhat in the optimal range
	
		// params
		const float MAX_GROW_FACTOR = 0.33f;				// maximum gain in total list length allowed by a node split operation

	private:
		struct coordPred
		{
			const SphereIndex* mSphereIndex;
			uint mDim;
			
			coordPred(const SphereIndex* sphereIndex, uint dim) : mSphereIndex(sphereIndex), mDim(dim) {}
			bool operator()(const uint& a, const uint& b) const { return mSphereIndex->getSphereCenterDim(a, mDim) < mSphereIndex->getSphereCenterDim(b, mDim); }
		};


		struct Node
		{
			Node() {}
			Node(vector<uint>* data_) : data(data_) {}
			Node(vector<uint>* data_, const BBox& bbox_) : data(data_), bbox(bbox_) {}
			~Node()
			{
				if (left)	delete left;
				if (right)	delete right;
				if (data)	delete[] data;
			}

			Node* left = nullptr;
			Node* right = nullptr;
			vector<uint>* data = nullptr;
			BBox bbox;
		};


	private:
		Node* mRoot = nullptr;

		const vector<vec4>* mSpheresVarying = nullptr;		// pointer for spheres with individual radius in w-element

		const vector<vec3>* mSpheresConstant = nullptr;		// pointer for spheres with constant radius stored in mSpheresConstantRadius
		float mSpheresConstantRadius = -1;
		bool mUsesConstRadius = false;
			
		
	public:
		SphereIndex()
		{
		}

		SphereIndex(const vector<vec4>& spheres, uint leafSize = DEFAULT_LEAF_SIZE)
		{
			create(spheres, leafSize);
		}

		SphereIndex(const vector<vec3>& sphereCenters, float radius, uint leafSize = DEFAULT_LEAF_SIZE)
		{
			create(sphereCenters, radius, leafSize);
		}

		~SphereIndex()
		{
			if (mRoot)
				delete mRoot;
		}


		// creates a kd-tree index on the given set of spheres.
		void create(const vector<vec3>& sphereCenters, float radius, uint leafSize = DEFAULT_LEAF_SIZE)
		{
			assert(!sphereCenters.empty());	// can't create index on empty set

			// release current index
			if (mRoot)
				delete mRoot;

			// save pointer to spheres data array
			mSpheresVarying = nullptr;
			mSpheresConstant = &sphereCenters;
			mSpheresConstantRadius = radius;
			mUsesConstRadius = true;

			buildUpKdTree(leafSize);
		}


		// creates a kd-tree index on the given set of spheres.
		void create(const vector<vec4>& spheres, uint leafSize = DEFAULT_LEAF_SIZE)
		{
			assert(!spheres.empty());	// can't create index on empty set

			// release current index
			if (mRoot)
				delete mRoot;

			// save pointer to spheres data array
			mSpheresConstant = nullptr;
			mSpheresConstantRadius = -1;
			mSpheresVarying = &spheres;
			mUsesConstRadius = false;

			buildUpKdTree(leafSize);
		}



		// queries all spheres in the index intersecting a given search ball.
		// all previous content in outIndices will be cleared.
		void radiusSearch(const vec3& queryPoint, float radius, vector<uint>& outIndices) const
		{
			outIndices.clear();

			// flag array indicating if a sphere was already found and added to outIndices
			vector<bool> found(mUsesConstRadius ? mSpheresConstant->size() : mSpheresVarying->size(), 0);

			// traverse tree
			vector<Node*> workingNodes = { mRoot };
			while (!workingNodes.empty())
			{
				// pick node from stack
				Node* node = workingNodes.back();
				workingNodes.pop_back();

				// test for bbox intersection
				if (node->bbox.intersectsSphere(queryPoint, radius))
				{
					// leaf node -> test individual spheres for intersection
					if (node->data)
					{
						for (uint sidx : *node->data)
						{
							if (found[sidx])
								continue;

							vec4 s = getSphere(sidx);
							if (Sphere(queryPoint, radius).intersects(Sphere(s.xyz(), s.w)))
							{
								outIndices.push_back(sidx);
								found[sidx] = 1;
							}
						}
					}
					// internal node -> push children on stack
					else
					{
						workingNodes.push_back(node->left);
						workingNodes.push_back(node->right);
					}
				}
			}
		}


		struct NeighborProcessor
		{
			virtual void operator() (uint nIndex, const vec3& npos) {}

			virtual void finalize() {}
		};


		void processNeighbors(const vec3& queryPoint, float radius, NeighborProcessor& nproc) const
		{
			// flag array indicating if a sphere was already found and added to outIndices
			vector<bool> found(mUsesConstRadius ? mSpheresConstant->size() : mSpheresVarying->size(), 0);

			// traverse tree
			vector<Node*> workingNodes = { mRoot };
			while (!workingNodes.empty())
			{
				// pick node from stack
				Node* node = workingNodes.back();
				workingNodes.pop_back();

				// test for bbox intersection
				if (node->bbox.intersectsSphere(queryPoint, radius))
				{
					// leaf node -> test individual spheres for intersection
					if (node->data)
					{
						for (uint sidx : *node->data)
						{
							if (found[sidx])
								continue;

							vec4 s = getSphere(sidx);
							if (Sphere(queryPoint, radius).intersects(Sphere(s.xyz(), s.w)))
							{
								nproc(sidx, s.xyz());
								found[sidx] = 1;
							}
						}
					}
					// internal node -> push children on stack
					else
					{
						workingNodes.push_back(node->left);
						workingNodes.push_back(node->right);
					}
				}
			}

			nproc.finalize();
		}


		void processNeighborsParallel(const vec3& queryPoint, float radius, NeighborProcessor& nproc) const
		{
			// flag array indicating in-range points
			vector<int> inside(mUsesConstRadius ? mSpheresConstant->size() : mSpheresVarying->size(), 0);
			
			// 1. traverse tree and mark all points inside the query radius
			vector<Node*> workingNodes = { mRoot };
			while (!workingNodes.empty())
			{
				// pick node from stack
				Node* node = workingNodes.back();
				workingNodes.pop_back();

				// test for bbox intersection
				if (node->bbox.intersectsSphere(queryPoint, radius))
				{
					// leaf node -> test individual spheres for intersection
					if (node->data)
					{
						#pragma omp parallel for
						for (int i = 0; i < (int)node->data->size(); i++)
						{
							uint sidx = node->data->at(i);
							vec4 s = getSphere(sidx);
							if (Sphere(queryPoint, radius).intersects(Sphere(s.xyz(), s.w)))
								inside[sidx] = 1;
						}
					}
					// internal node -> push children on stack
					else
					{
						workingNodes.push_back(node->left);
						workingNodes.push_back(node->right);
					}
				}
			}


			// 2. process selected neighbors in parallel
			for (int i = 0; i < (int)inside.size(); i++)
				if (inside[i])
					nproc(i, getSphere(i).xyz());

			nproc.finalize();
		}



		// returns a vector of pointers to the index arrays of the kd-tree's leaf-nodes
		void getLeafIndices(vector<const vector<uint>*>& outLeafIndices) const
		{
			outLeafIndices.clear();

			// traverse tree
			vector<Node*> workingNodes = { mRoot };
			while (!workingNodes.empty())
			{
				// pick node from stack
				Node* node = workingNodes.back();
				workingNodes.pop_back();

				if (node->data)
				{
					// leaf node
					outLeafIndices.push_back(node->data);
				}
				else
				{
					// internal node -> push children on stack
					workingNodes.push_back(node->left);
					workingNodes.push_back(node->right);
				}
			}
		}


	protected:
		// builds up the kd-tree after the sphere data was set to the sphere arrays, independent of the internal representation.
		void buildUpKdTree(uint leafSize = DEFAULT_LEAF_SIZE)
		{
			uint size = uint(mUsesConstRadius ? mSpheresConstant->size() : mSpheresVarying->size());

			// compute bounding box
			BBox bbox;
			for (uint i = 0; i < size; i++)
			{
				const vec4& s = getSphere(i);
				bbox.extend(BBox(s.xyz() - vec3(s.w, s.w, s.w), s.xyz() + vec3(s.w, s.w, s.w)));
			}

			// create index list and sort in x,y,z direction by sphere centers
			vector<uint>* lists = new vector<uint>[3];
			lists[0].resize(size);
			for (uint i = 0; i < size; ++i)	lists->at(i) = i;	// fill list 0
			lists[2] = lists[1] = lists[0];						// copy to list 1 and 2
			// sort
			for (uint dim = 0; dim < 3; ++dim)
				parallel::sort(lists[dim], coordPred(this, dim));	// faster than serial sort
			

			// insert root node
			mRoot = new Node(lists, bbox);

			// recursively split
			vector<Node*> workingNodes = { mRoot };
			while (!workingNodes.empty())
			{
				// pick next working node
				Node* node = workingNodes.back();
				workingNodes.pop_back();

				// too large for leaf node -> split
				if (node->data->size() > leafSize)
				{
					// determine longest bb side
					vec3 sides = fabsf(node->bbox.dim());
					int dim = sides.x > sides.y ? (sides.x > sides.z ? 0 : 2) : (sides.y > sides.z ? 1 : 2);

					// determine splitting plane in dim. The volume is always split by one of its sphere centers
					uint splitIndex = node->data->size() / 2;
					uint splitSidx = node->data[dim][splitIndex];
					float splitCoord = getSphereCenterDim(splitSidx, dim);

					// split the sorted lists by the relative location of the sphere centers to the splitting plane
					Node* leftChild = new Node(new vector<uint>[3]);
					Node* rightChild = new Node(new vector<uint>[3]);
					for (uint d = 0; d < 3; ++d)
					{
						const vector<uint> &dimList = node->data[d];

						// split list along dimension d
						for (uint sidx : dimList)
						{
							const vec4& s = getSphere(sidx);
							BBox sbbox(s.xyz() - vec3(s.w, s.w, s.w), s.xyz() + vec3(s.w, s.w, s.w));

							// test sphere intersection left volume / right volume
							// update left/right bounding box only once  
							if (s[dim] - s.w < splitCoord)
							{
								leftChild->data[d].push_back(sidx);
								if (d == 0) leftChild->bbox.extend(sbbox);
							}
							if (s[dim] + s.w >= splitCoord)
							{
								rightChild->data[d].push_back(sidx);
								if (d == 0) rightChild->bbox.extend(sbbox);
							}
						}
					}


					float growFactor = float(leftChild->data->size() + rightChild->data->size()) / float(node->data->size()) - 1.0f;
					if (growFactor > MAX_GROW_FACTOR)
					{
						//cout << "Node split increase " << (leftChild->data->size() + rightChild->data->size()) << "/" << node->data->size() << "   (" << (growFactor * 100) << "%)" << endl;
					}


					// if the current split has reduced the list sizes, proceed
					if (leftChild->data->size() < node->data->size() && rightChild->data->size() < node->data->size() && growFactor <= MAX_GROW_FACTOR)
					{
						// delete data list, adopt children and push them on stack
						node->left = leftChild;
						node->right = rightChild;
						delete[] node->data;
						node->data = nullptr;

						workingNodes.push_back(node->left);
						workingNodes.push_back(node->right);
					}
					else
					{
						// node remains leaf node
						delete leftChild;
						delete rightChild;
						node->data[1].clear();
						node->data[2].clear();
					}
				}
				// leaf node
				else
				{
					// saving indices 3 times not necessary - just keep x-direction
					node->data[1].clear();
					node->data[2].clear();
				}

			} // end while 
		}



		// returns the sphere with the given index, independent of the internal representation (constant or varying radii)
		vec4 getSphere(uint index) const
		{
			if (mUsesConstRadius)
				return vec4(mSpheresConstant->at(index), mSpheresConstantRadius);
			else
				return mSpheresVarying->at(index);
		}

		// returns the dimension dim of the sphere center with the given index, independent of the internal representation (constant or varying radii)
		const float getSphereCenterDim(uint index, uint dim) const
		{
			if (mUsesConstRadius)
				return mSpheresConstant->at(index)[dim];
			else
				return mSpheresVarying->at(index)[dim];
		}

	};	/// class SphereIndex


}	/// end namespace gms

