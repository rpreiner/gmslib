//-----------------------------------------------------------------------------
// gmslib - Gaussian Mixture Surface Library
// Copyright (c) Reinhold Preiner 2014-2020 
// 
// Usage is subject to the terms of the WFP (modified BSD-3-Clause) license.
// See the accompanied LICENSE file or
// https://github.com/rpreiner/gmslib/blob/main/LICENSE
//-----------------------------------------------------------------------------

#pragma once

#pragma warning(disable: 4996)
#include "ext/nanoply.hpp"
#include "vec.hpp"
#include "pointset.hpp"
#include "mixture.hpp"


namespace gms
{

	class IO
	{
	public:
		static bool loadPointCloudPLY(const string& filename, gms::PointSet* points)
		{
			//Get file info
			nanoply::Info info(filename.c_str());
			if (info.errInfo != nanoply::NNP_OK)
			{
				std::cout << "Invalid file format" << std::endl;
				return false;
			}

			//Resize the element containers
			int vertCnt = info.GetVertexCount();
			if (vertCnt <= 0)
			{
				std::cout << "The file does't contain any vertex." << std::endl;
				return false;
			}

			// prepare buffers
			points->resize(vertCnt);
			shared_ptr<vector<gms::vec3>> normalArray = nullptr;
			shared_ptr<vector<gms::vec4>> ellipseArray = nullptr;
			shared_ptr<vector<float>> qualityArray = nullptr;
			shared_ptr<vector<gms::vec3>> colorArray = nullptr;

			//Create the vertex properties descriptor (what ply property and where to save its data)
			nanoply::ElementDescriptor vertex(nanoply::NNP_VERTEX_ELEM);
			if (vertCnt > 0)
			{
				vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<gms::vec3, 3, float>(nanoply::NNP_PXYZ, points->data()));

				// determine the color format of this ply file
				nanoply::PlyEntity colorEntity = nanoply::NNP_CRGBA;
				for (const auto& elem : info.elemVec)
					for (const auto& prop : elem.propVec)
					{
						if (prop.name == "rgb")
						{
							colorEntity = prop.elem;
							colorArray = make_shared<vector<gms::vec3>>(points->size());
							vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<gms::vec3, 3, float>(colorEntity, colorArray->data()));
						}
						else if (prop.elem == nanoply::NNP_NXYZ)
						{
							normalArray = make_shared<vector<gms::vec3>>(points->size());
							vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<gms::vec3, 3, float>(nanoply::NNP_NXYZ, normalArray->data()));
						}
						else if (prop.elem == nanoply::NNP_QUALITY)
						{
							qualityArray = make_shared<vector<float>>(points->size());
							vertex.dataDescriptor.push_back(new nanoply::DataDescriptor<float, 1, float>(nanoply::NNP_QUALITY, qualityArray->data()));
						}
					}
			}

			//Create the mesh descriptor
			std::vector<nanoply::ElementDescriptor*> meshDescr;
			meshDescr.push_back(&vertex);

			//Open the file and save the element data according the relative element descriptor
			OpenModel(info, meshDescr);
			for (int i = 0; i < vertex.dataDescriptor.size(); i++)
				delete vertex.dataDescriptor[i];

			if (normalArray != nullptr)		points->setAttrib("normal", normalArray);
			if (qualityArray != nullptr)	points->setAttrib("quality", qualityArray);
			if (colorArray != nullptr)		points->setAttrib("color", colorArray);

			return (info.errInfo == nanoply::NNP_OK);
		}


		static bool saveMixturePLY(const string& filename, const gms::Mixture& mixture)
		{
			string plyString;
			exportMixture(mixture, plyString);

			// write to file
			std::ofstream out(filename);
			out << plyString;
			out.close();

			if (!out)
			{
				cerr << "Error writing to file " << filename << endl;
				return false;
			}
			cout << "Wrote mixture to file " << filename << endl;
			return true;
		}


		/// export mixture to ascii PLY file format string
		static void exportMixture(const gms::Mixture& mixture, string& plyString)
		{
			stringstream ply;
			ply << "ply" << endl;
			ply << "format ascii 1.0" << endl;
			ply << "comment covmesh 2.0" << endl;
			ply << "comment gms exported Gaussian mixture ply file" << endl;
			ply << "comment (c) Reinhold Preiner" << endl;

			ply << "element component " << mixture.size() << endl;
			ply << "property float x" << endl;
			ply << "property float y" << endl;
			ply << "property float z" << endl;
			ply << "property float covxx" << endl;
			ply << "property float covxy" << endl;
			ply << "property float covxz" << endl;
			ply << "property float covyy" << endl;
			ply << "property float covyz" << endl;
			ply << "property float covzz" << endl;
			ply << "property float weight" << endl;
			if (mixture.hasNormals())
			{
				ply << "property float nvx" << endl;
				ply << "property float nvy" << endl;
				ply << "property float nvz" << endl;
			}
			ply << "end_header" << endl;

			// write components
			for (unsigned i = 0; i < mixture.size(); ++i)
			{
				const gms::Gaussian& g = mixture[i];
				ply << g.mu.x << "  ";
				ply << g.mu.y << "  ";
				ply << g.mu.z << "  ";
				ply << g.cov.e00 << "  ";
				ply << g.cov.e01 << "  ";
				ply << g.cov.e02 << "  ";
				ply << g.cov.e11 << "  ";
				ply << g.cov.e12 << "  ";
				ply << g.cov.e22 << "  ";
				ply << g.weight << "  ";
				if (mixture.hasNormals())
				{
					auto& nvar = mixture.nvars[i];
					ply << nvar.x << "  ";
					ply << nvar.y << "  ";
					ply << nvar.z << "  ";
				}
				ply << endl;
			}

			plyString = ply.str();
		}


		/// load mixture from ascii PLY file format string.
		/// if non-empty mixture is passed, the loaded components are directly appended to the existing ones.
		static void importMixture(const string& plyString, gms::Mixture& mixture)
		{
			typedef enum
			{
				X, Y, Z, COVXX, COVXY, COVXZ, COVYY, COVYZ, COVZZ, W, NVX, NVY, NVZ

			} MixtureElem;

			string line;
			int componentCounter = 0;	// running counter
			int componentCount = 0;		// number of components
			bool readData = false;
			int lineCount = 0;
			bool hasNormals = false;
			vector<MixtureElem> header;

			// breacket first line
			string::size_type lineStart = 0;
			string::size_type lineEnd = plyString.find_first_of("\n", lineStart);

			while (lineStart != string::npos)
			{
				// extract line
				line = plyString.substr(lineStart, lineEnd - lineStart);

				// tokenize line
				string::size_type pos1 = line.find_first_not_of(" \n\r\t", 0);
				string::size_type pos2 = line.find_first_of(" \n\r\t", pos1);
				if (pos1 != string::npos && line.substr(pos1, pos2 - pos1).compare("comment") != 0)
				{
					// "element" command			
					if (line.substr(pos1, pos2 - pos1).compare("element") == 0)
					{
						pos1 = line.find_first_not_of(" \n\r\t", pos2);
						pos2 = line.find_first_of(" \n\r\t", pos1);

						// Gaussian components
						if (line.substr(pos1, pos2 - pos1).compare("component") == 0)
						{
							pos1 = line.find_first_not_of(" \n\r\t", pos2);
							pos2 = line.find_first_of(" \n\r\t", pos1);
							componentCount = atoi(line.substr(pos1, pos2 - pos1).c_str());
						}
					}
					// "property" command			
					else if (line.substr(pos1, pos2 - pos1).compare("property") == 0)
					{
						// read "property type"
						pos1 = line.find_first_not_of(" \n\r\t", pos2);
						pos2 = line.find_first_of(" \n\r\t", pos1);

						// don't need to parse the list property
						if (line.substr(pos1, pos2 - pos1).compare("list") != 0)
						{
							// read property name
							pos1 = line.find_first_not_of(" \n\r\t", pos2);
							pos2 = line.find_first_of(" \n\r\t", pos1);
							string propName = line.substr(pos1, pos2 - pos1);

							// save all elements that are contained in this file in header vector
							MixtureElem elem;
							if (propName.compare("x") == 0) elem = X;
							else if (propName.compare("y") == 0) elem = Y;
							else if (propName.compare("z") == 0) elem = Z;
							else if (propName.compare("covxx") == 0) elem = COVXX;
							else if (propName.compare("covxy") == 0) elem = COVXY;
							else if (propName.compare("covxz") == 0) elem = COVXZ;
							else if (propName.compare("covyy") == 0) elem = COVYY;
							else if (propName.compare("covyz") == 0) elem = COVYZ;
							else if (propName.compare("covzz") == 0) elem = COVZZ;
							else if (propName.compare("weight") == 0) elem = W;
							else if (propName.compare("nvx") == 0) { elem = NVX; hasNormals = true; }
							else if (propName.compare("nvy") == 0) elem = NVY;
							else if (propName.compare("nvz") == 0) elem = NVZ;
							header.push_back(elem);
						}
					}
					// "end_header"
					else if (line.substr(pos1, pos2 - pos1).compare("end_header") == 0)
					{
						readData = true;
					}
					// element data
					else if (readData)
					{
						// parse components
						if (componentCounter < componentCount)
						{
							gms::vec3 mu, nvar;
							gms::smat3 cov;
							float w;

							for (int i = 0; pos1 != string::npos; ++i)
							{
								float val = (float)atof(line.substr(pos1, pos2 - pos1).c_str());

								if (i < header.size()) switch (header[i])
								{
								case X: mu.x = val; break;
								case Y: mu.y = val; break;
								case Z: mu.z = val; break;
								case COVXX: cov.e00 = val; break;
								case COVXY: cov.e01 = val; break;
								case COVXZ: cov.e02 = val; break;
								case COVYY: cov.e11 = val; break;
								case COVYZ: cov.e12 = val; break;
								case COVZZ: cov.e22 = val; break;
								case W: w = val; break;
								case NVX: nvar.x = val; break;
								case NVY: nvar.y = val; break;
								case NVZ: nvar.z = val; break;
								}

								pos1 = line.find_first_not_of(" \n\r\t", pos2);
								pos2 = line.find_first_of(" \n\r\t", pos1);
							}

							mixture.push_back(gms::Gaussian(mu, cov, w));

							if (hasNormals)	mixture.nvars.push_back(nvar);

							++componentCounter;
						}
					}
				}

				// bracket next line
				lineStart = plyString.find_first_not_of("\n", lineEnd);
				lineEnd = plyString.find_first_of("\n", lineStart);
			}
		}

	};

} // end namespace gms
