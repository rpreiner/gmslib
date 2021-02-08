//-----------------------------------------------------------------------------
// gmslib - Gaussian Mixture Surface Library
// (c) Reinhold Preiner 2014-2020 
// 
// Usage is subject to the terms of the WFP (modified BSD-3-Clause) license.
// See the accompanied LICENSE file or
// https://github.com/rpreiner/gmslib/blob/main/LICENSE
//-----------------------------------------------------------------------------

#include "gmslib/argparser.hpp"
#include "gmslib/pointset.hpp"
#include "gmslib/mixture.hpp"
#include "gmslib/io.hpp"
#include <iostream>
using namespace std;
using namespace gms;



int main(int argc, char* argv[])
{
	ArgParser parser;
	parser.addArgument("i", "Input point cloud file");
	parser.addArgument("o", "Output mixture file");

	parser.addArgument("verbose", "Verbose output");
	parser.addArgument("memory", "Activate memory usage profiling");
	parser.addArgument("alpha", "Clustering regularization parameter");
	parser.addArgument("blocksize", "Compute mixture in blocks of the specified point count");
	parser.addArgument("pointpos", "Initializes Gaussian positions in point locations rather than local point means");
	parser.addArgument("stdev", "Default isotropic standard deviation bias of each initial Gaussian [in %bbd]");
	parser.addArgument("iso", "Initialize mixture with isotropic Gaussians of standard deviation <stdev>");
	parser.addArgument("inittype", "'knn' - Init anisotropic Gaussians based on KNN; 'fixed' - based on fixed distance");
	parser.addArgument("knn", "Number of nearest neighbors considered for 'knn' initialization");
	parser.addArgument("fixeddist", "Max neighborhood distance for points considered for 'fixed' initialization [in %bbd]");
	parser.addArgument("weighted", "Initializes mixture with locally normalized density");
	parser.addArgument("levels", "Number of HEM clustering levels");
	parser.addArgument("threads", "Number of parallel threads");
	parser.addNote("Quantities described with '[in %bbd]' are given in percent of the input point cloud bounding box diagonal.");

	if (argc == 1 || (argc > 1 && string(argv[1]) == "-help"))
	{
		cout << parser.helpStr() << endl;
		return 0;
	}
	if (!parser.matchArgs(argc, argv))
		return 1;
	//------------------------------------------------------------------

	string infilename = parser.getArgument("i");
	string outfilename = parser.getArgument("o");
	if (infilename == "") { cerr << "No input point cloud file specified." << endl; return 1; }
	if (outfilename == "") { cerr << "No output point cloud file specified." << endl; return 1; }


	Mixture::Params params;
	params.computeNVar = false;		// deactivate CLOP normal clustering 
	params.verbose = parser.getBool("verbose", true);
	params.memoryProfiling = parser.getBool("memory", false);
	params.alpha = parser.getFloat("alpha", 2.0f);
	params.blockSize = parser.getUint("blocksize", 0);
	params.blockProcessing = params.blockSize > 0;
	params.hemReductionFactor = 3.0f;
	params.initIsotropic = parser.getBool("iso", false);
	params.initIsotropicStdev = parser.getFloat("stdev", 0.01f);
	params.initMeansInPoints = parser.getBool("pointpos", true);
	params.kNNCount = parser.getUint("knn", 8);
	params.maxInitNeighborDist = parser.getFloat("fixeddist", 0.1f);
	params.nLevels = parser.getUint("levels", 20);
	params.numThreads = parser.getUint("threads", 8);
	params.useWeightedPotentials = parser.getBool("weighted", false);
	params.initNeighborhoodType = 0;
	string s = parser.getArgument("inittype");
	if (s != "")
	{
		if (s == "fixed")		params.initNeighborhoodType = 0;
		else if (s == "knn")	params.initNeighborhoodType = 1;
		else {
			std::cerr << "Invalid 'anisotype' argument. Use 'knn' or 'fixed'.\n";
			exit(1);
		}
	}
	//-----------------------------------------------------------------------

	auto points = new gms::PointSet();
	if (IO::loadPointCloudPLY(infilename, points) == false) {
		std::cerr << "Cannot load point cloud file '" << infilename << "'\n";
		exit(1);
	}


	// convert bbd-relative quantities
	gms::BBox bboxPoints(*points);
	float conversionFac = 0.01f * length(bboxPoints.dim());
	params.maxInitNeighborDist *= conversionFac;
	params.initIsotropicStdev *= conversionFac;

	// compute mixture
	cout << "\n--- Compute Mixture ---\n";
	Timer timer;
	Mixture mixture(points, params);
	double t = timer.stop() * 0.001;

	cout << endl << std::fixed << std::setprecision(5) << "Created mixture in " << t << " s" << endl;
	cout << "# Gaussians: " << mixture.size() << endl << endl;

	// write output
	if (!IO::saveMixturePLY(outfilename, mixture))
		exit(1);

	return 0;
}
