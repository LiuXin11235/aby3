
#include <cryptoTools/Common/CLP.h>
#include <map>

#include "aby3_tests/aby3_tests.h"
#include "eric.h"
#include "aby3-DB-main.h"
#include "aby3-DB_tests/UnitTests.h"
#include <tests_cryptoTools/UnitTests.h>
#include <aby3-ML/main-linear.h>
#include <aby3-ML/main-logistic.h>
#include <aby3-RTR/CipherIndex.h>

#include "tests_cryptoTools/UnitTests.h"
#include "cryptoTools/Crypto/PRNG.h"

using namespace oc;
using namespace aby3;
std::vector<std::string> unitTestTag{ "u", "unitTest" };

#define BASIC_TEST

void help()
{

	std::cout << "-u                        ~~ to run all tests" << std::endl;
	std::cout << "-u n1 [n2 ...]            ~~ to run test n1, n2, ..." << std::endl;
	std::cout << "-u -list                  ~~ to list all tests" << std::endl;
	std::cout << "-intersect -nn NN [-c C]  ~~ to run the intersection benchmark with 2^NN set sizes, C 32-bit data columns." << std::endl;
	std::cout << "-eric -nn NN              ~~ to run the eric benchmark with 2^NN set sizes" << std::endl;
	std::cout << "-threat -nn NN -s S       ~~ to run the threat log benchmark with 2^NN set sizes and S sets" << std::endl;
}


// int main(int argc, char** argv)
// {


// 	try {


// 		bool set = false;
// 		oc::CLP cmd(argc, argv);


// 		if (cmd.isSet(unitTestTag))
// 		{
// 			auto tests = tests_cryptoTools::Tests;
// 			tests += aby3_tests;
// 			tests += DB_tests;

// 			tests.runIf(cmd);
// 			return 0;
// 		}

// 		if (cmd.isSet("RTR"))
// 		{
// 			set = true;
// 			test_mul(cmd);
// 			std::cout << "can call functions " << std::endl;
// 		}

// 		if (cmd.isSet("linear-plain"))
// 		{
// 			set = true;
// 			linear_plain_main(cmd);
// 		}
// 		if (cmd.isSet("linear"))
// 		{
// 			set = true;
// 			linear_main_3pc_sh(cmd);
// 		}

// 		if (cmd.isSet("logistic-plain"))
// 		{
// 			set = true;
// 			logistic_plain_main(cmd);
// 		}

// 		if (cmd.isSet("logistic"))
// 		{
// 			set = true;
// 			logistic_main_3pc_sh(cmd);
// 		}

// 		if (cmd.isSet("eric"))
// 		{
// 			set = true;

// 			auto nn = cmd.getMany<int>("nn");
// 			if (nn.size() == 0)
// 				nn.push_back(16);

// 			for (auto n : nn)
// 			{
// 				eric(1 << n);
// 			}
// 		}


// 		if (cmd.isSet("intersect"))
// 		{
// 			set = true;

// 			auto nn = cmd.getMany<int>("nn");
// 			auto c = cmd.getOr("c", 0);
// 			if (nn.size() == 0)
// 				nn.push_back(1 << 16);

// 			for (auto n : nn)
// 			{
// 				auto size = 1 << n;
// 				DB_Intersect(size, c, cmd.isSet("sum"));
// 			}
// 		}


// 		if (cmd.isSet("threat"))
// 		{
// 			set = true;

// 			auto nn = cmd.getMany<int>("nn");
// 			auto c = cmd.getOr("s", 2);
// 			if (nn.size() == 0)
// 				nn.push_back(1 << 16);

// 			for (auto n : nn)
// 			{
// 				auto size = 1 << n;
// 				DB_threat(size, c);
// 			}
// 		}



// 		if (cmd.isSet("card"))
// 		{
// 			set = true;

// 			auto nn = cmd.getMany<int>("nn");
// 			if (nn.size() == 0)
// 				nn.push_back(1 << 16);

// 			for (auto n : nn)
// 			{
// 				auto size = 1 << n;
// 				DB_cardinality(size);
// 			}
// 		}
		
// 		//if (cmd.isSet("add"))
// 		//{
// 		//	set = true;

// 		//	auto nn = cmd.getMany<int>("nn");
// 		//	if (nn.size() == 0)
// 		//		nn.push_back(1 << 16);

// 		//	for (auto n : nn)
// 		//	{
// 		//		auto size = 1 << n;
// 		//		Sh3_add_test(size);
// 		//	}
// 		//}

// 		if (set == false)
// 		{
// 			help();
// 		}

// 	}
// 	catch (std::exception& e)
// 	{
// 		std::cout << e.what() << std::endl;
// 	}

// 	return 0;
// }

int main(int argc, char** argv)
{
	oc::CLP cmd(argc, argv);

	#ifdef PERFORMANCE_TEST
	// test the vectorization for basic ops (mul) and (gt).
	int repeats = 20;
	int step = 50, start = 1, end = 1e3;
	int points = (end - start) / step;
	std::vector<int> n_list;
	for(int i=0; i<points; i++){
		n_list.push_back(start + i*step);
	}

	std::map<int, std::map<std::string, std::vector<double>>> performance_dict;
	for(int i=0; i<n_list.size(); i++){
		int n = n_list[i];
		std::map<std::string, std::vector<double>> tmp_map;
		basic_performance(cmd, n, repeats, tmp_map);
		performance_dict[n] = tmp_map;
	}

	// cout the result.
	for(int i=0; i<n_list.size(); i++){
		std::cout << "\nvector size = " << n_list[i] << std::endl;
		std::cout << "mul" << std::endl;
		for(int j=0; j<3; j++){
			std::cout << performance_dict[n_list[i]]["mul"][j] << " ";
		}
		std::cout << "\ngt" << std::endl;
		for(int j=0; j<3; j++){
			std::cout << performance_dict[n_list[i]]["gt"][j] << " ";
		}
		std::cout << std::endl;
	}
	#endif

	#ifdef BASIC_TEST
	// test gt
	test_gt(cmd);
	// test eq
	// test_eq(cmd);

	// test cipher_argsort
	test_argsort(cmd);
	#endif

	return 0;
}
