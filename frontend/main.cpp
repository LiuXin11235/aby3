
#include <cryptoTools/Common/CLP.h>
#include <map>

#include "eric.h"
#include <aby3-ML/main-linear.h>
#include <aby3-ML/main-logistic.h>
#include "aby3-RTR/RTRTest.h"
#include "aby3-RTR/DistributeRTRTest.h"


using namespace oc;
using namespace aby3;
// std::vector<std::string> unitTestTag{ "u", "unitTest" };

int main(int argc, char** argv)
{
	oc::CLP cmd(argc, argv);

	int prog = -1;
	if(cmd.isSet("prog")){
		auto progs = cmd.getMany<int>("prog");
		prog = progs[0];
	}
	if(prog == -1){
		std::cout << "Set prog to 0(basic performance test) by default " << std::endl;
	}

	if(prog == 0){ // test the vectorization for basic ops (mul) and (gt).
		int repeats = int(100);

		std::vector<int> n_list = {      10,       13,       17,       23,       30,       40,
				54,       71,       95,      126,      167,      222,
				294,      390,      517,      686,      910,     1206,
			1599,     2120,     2811,     3727,     4941,     6551,
			8685,    11513,    15264,    20235,    26826,    35564,
			47148,    62505,    82864,   109854,   145634,   193069,
			255954,   339322,   449843,   596362,   790604,  1048113,
			1389495,  1842069,  2442053,  3237457,  4291934,  5689866,
			7543120, 10000000};

		std::map<int, std::map<std::string, std::vector<double>>> performance_dict;
		for(int i=0; i<n_list.size(); i++){
			int n = n_list[i];
			std::map<std::string, std::vector<double>> tmp_map;
			basic_performance(cmd, n, repeats, tmp_map);
			// dis_basic_performance(cmd, n, repeats, tmp_map);
			performance_dict[n] = tmp_map;

			// execute one evaluation and record one.
			std::cout << "\nvector size = " << n_list[i] << std::endl;
			std::cout << "mul" << std::endl;
			for(int j=0; j<3; j++){
				std::cout << performance_dict[n_list[i]]["mul"][j] << " ";
			}
			std::cout << "\ngt" << std::endl;
			for(int j=0; j<3; j++){
				std::cout << performance_dict[n_list[i]]["gt"][j] << " ";
			}
			// std::cout << std::endl;
			std::cout << "\nadd" << std::endl;
			for(int j=0; j<3; j++){
				std::cout << performance_dict[n_list[i]]["add"][j] << " ";
			}
			std::cout << std::endl;
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
			// std::cout << std::endl;
			std::cout << "\nadd" << std::endl;
			for(int j=0; j<3; j++){
				std::cout << performance_dict[n_list[i]]["add"][j] << " ";
			}
			std::cout << std::endl;
		}

		return 0;
	}

	if(prog == 1){ // test the performance in the distribued setting.
		int repeats;
		std::vector<int> n_list = {      10,       13,       17,       23,       30,       40,
					54,       71,       95,      126,      167,      222,
				294,      390,      517,      686,      910,     1206,
				1599,     2120,     2811,     3727,     4941,     6551,
				8685,    11513,    15264,    20235,    26826,    35564,
				47148,    62505,    82864,   109854,   145634,   193069,
				255954,   339322,   449843,   596362,   790604,  1048113,
			1389495,  1842069,  2442053,  3237457,  4291934,  5689866,
			7543120, 10000000};

		std::map<int, std::map<std::string, double>> performance_dict;
		for(int i=0; i<n_list.size(); i++){

			int n = n_list[i];
			
			// set the repeat times.
			if(i < 20) repeats = int(1e4);
			else if(i < 40) repeats = int(1e3);
			else repeats = int(100);

			std::map<std::string, double> tmp_map;
			dis_basic_performance(cmd, n, repeats, tmp_map);
			performance_dict[n] = tmp_map;

			// execute one evaluation and record one.
			std::cout << "\nvector size = " << n_list[i] << std::endl;
			std::map<std::string, double>::iterator iter;
			iter = tmp_map.begin();
			while(iter != tmp_map.end()){
				std::cout << iter->first << std::endl;
				std::cout << iter->second << std::endl;
				iter ++;
			}
		}
		return 0;
	}

	if(prog == 2){ // basic test.
		// test gt
		test_gt(cmd);

		// test eq - has problems.
		test_eq(cmd);

		// test multiplication between bits and ints.
		test_mul(cmd);

		// test cipher_argsort
		test_argsort(cmd, 1);
		test_argsort(cmd, 0);

		// test cipher_index
		test_cipher_index(cmd, 0);
		test_cipher_index(cmd, 1);

		// test binning.
		test_cipher_binning(cmd, 0);
		test_cipher_binning(cmd, 1);

		return 0;
	}

	if(prog == 3){ // basic test in the distributed setting.
		// test_mul(cmd);
		// dis_test_mul(cmd);
		int repeats = int(10);

		std::vector<int> n_list = {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960, 81920, 163840, 327680, 655360, 1310720, 2621440, 5242880, 10485760, 20971520};
		std::vector<int> m_list;

		if(cmd.isSet("mRatio")){
			auto m_ratios = cmd.getMany<double>("mRatio");
			double m_ratio = m_ratios[0];
			for(int i=0; i<n_list.size(); i++){
				// m_list.push_back(int(m_ratio * n_list[i]));
				m_list.push_back(n_list[i] * m_ratio);
			}
		}
		else{
			int m_value = 10;
			if(cmd.isSet("mValue")){
				auto m_values = cmd.getMany<double>("mValue");
				m_value = m_values[0];
			}
			for(int i=0; i<n_list.size(); i++){
				// m_list.push_back(int(m_ratio * n_list[i]));
				m_list.push_back(m_value);
			}
		}

		int testFlag = -1;
		if(cmd.isSet("testFlag")){
			auto testFlags = cmd.getMany<int>("testFlag");
			testFlag = testFlags[0];
			if(testFlag > 2) {
				std::cout << "testFlag should be within 0-2" << std::endl;
				testFlag = -1;
			}
		}
		// test and output the time.
		for(int i=0; i<n_list.size(); i++){
			std::map<std::string, double> dict;
			dis_cipher_index_performance(cmd, n_list[i], m_list[i], repeats, dict, testFlag);
			std::cout << "n = " << n_list[i] << " m = " << m_list[i] << " normal_time = " << dict["normal"] << " plain_time = " << dict["plain"] << " rtr_time = " << dict["rtr"] << std::endl;
		}
		return 0;
	}
	
	if(prog == 4){ // test the functions in the new API.
		test_argsort(cmd, 2);
		test_cipher_index(cmd, 2);
	}

	// if(prog == 5){
	// 	test_cipher_index(cmd, 2);
	// }

	// if(prog == 6){
	// 	test_cipher_binning(cmd, 2);
	// }


	std::cout << "prog only support 0 - 4" << std::endl;
	return 0;
}
