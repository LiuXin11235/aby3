#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>

#include <time.h>

#include <cryptoTools/Network/IOService.h>

#include <cryptoTools/Common/CLP.h>
#include <Eigen/Dense>

#include "aby3-ML/Regression.h"
#include "aby3-ML/LinearModelGen.h"

#include "aby3-ML/aby3ML.h"
#include "aby3-ML/PlainML.h"


using namespace std;
using namespace Eigen;
using namespace oc;
namespace aby3
{

	bool readCSV(const std::string& filename, eMatrix<double>& matrix) {
			std::ifstream file(filename);
			std::cout << "Read " << filename << std::endl;
			if (!file.is_open()) {
					std::cerr << "Could not open the file: " << filename << std::endl;
					return false;
			}

			std::string line;
			std::vector<std::vector<double>> data;
			while (std::getline(file, line)) {
					std::stringstream lineStream(line);
					std::string cell;
					std::vector<double> row;
					while (std::getline(lineStream, cell, ',')) {
							row.push_back(std::stod(cell));
					}
					data.push_back(row);
			}

			if (data.empty()) {
					std::cerr << "No data found in the file: " << filename << std::endl;
					return false;
			}

			matrix.resize(data.size(), data[0].size());
			for (size_t i = 0; i < data.size(); ++i) {
					for (size_t j = 0; j < data[i].size(); ++j) {
							matrix(i, j) = data[i][j];
					}
			}

			return true;
	}
	
	// void printRevealedValues(const eMatrix<double>& W2, const LogisticModelGen& gen, bool print, u64 dim, aby3ML& p) {
	// 		auto w2Val = p.reveal(W2);

	// 		if (print) {
	// 				for (u64 i = 0; i < dim; ++i) {
	// 						std::cout << i << " " << gen.mModel(i, 0) << " " << w2Val(i, 0) << std::endl;
	// 				}
	// 		}
	// }

	int logistic_plain_main(CLP& cmd)
	{
		auto N = cmd.getOr<int>("N", 10000);
		auto D = cmd.getOr<int>("D", 1000);
		auto B = cmd.getOr<int>("B", 128);
		auto IT = cmd.getOr<int>("I", 10000);
		auto testN = cmd.getOr<int>("testN", 1000);

		PRNG prng(toBlock(1));
		LogisticModelGen gen;

		eMatrix<double> model(D, 1);
		for (u64 i = 0; i < (u64)D; ++i)
		{
			model(i, 0) = prng.get<int>() % 10;
			std::cout << model(i, 0) << " ";
		}

		std::cout << std::endl;
		gen.setModel(model);



		eMatrix<double> train_data(N, D), train_label(N, 1);
		eMatrix<double> test_data(testN, D), test_label(testN, 1);
		gen.sample(train_data, train_label);
		gen.sample(test_data, test_label);


		std::cout << "training __" << std::endl;

		RegressionParam params;
		params.mBatchSize = B;
		params.mIterations = IT;
		params.mLearningRate = 1.0 / (1 << 3);
		PlainML engine;

		eMatrix<double> W2(D, 1);
		W2.setZero();
		//engine.mPrint = cmd.isSet("print");

		SGD_Logistic(params, engine, train_data, train_label, W2, &test_data, &test_label);

		for (u64 i = 0; i < (u64)D; ++i)
		{
			std::cout << i << " " << gen.mModel(i, 0) << " " << W2(i, 0) << std::endl;
		}

		return 0;
	}


	int logistic_main_3pc_sh(int N, int dim, int B, int IT, int testN, int pIdx, bool print, CLP& cmd, Session& chlPrev, Session& chlNext, double thre = 0.5, std::string prefix = "lr_train/")
	{

		PRNG prng(toBlock(1));
		LogisticModelGen gen;

		eMatrix<double> val_train_data(N, dim), val_train_label(N, 1);
		eMatrix<double> val_test_data(testN, dim), val_test_label(testN, 1);
		// gen.sample(val_train_data, val_train_label);
		// gen.sample(val_test_data, val_test_label);

    // std::string prefix = "lr_train_toy/";
    if (!readCSV(prefix + "train_data.csv", val_train_data) ||
        !readCSV(prefix + "train_label.csv", val_train_label) ||
        !readCSV(prefix + "test_data.csv", val_test_data) ||
        !readCSV(prefix + "test_label.csv", val_test_label)) {
        std::cerr << "Error reading data from files." << std::endl;
        return -1;
    }		
    N = val_train_data.rows();
    dim = val_train_data.cols();
    testN = val_test_data.rows();

		eMatrix<double> val_W2(dim, 1);
		val_W2.setZero();

		eMatrix<double> model(dim, 1);
		for (u64 i = 0; i < (u64)std::min(dim, 10); ++i)
		{
			model(i, 0) = prng.get<int>() % 10;
		}
		gen.setModel(model);

		RegressionParam params;
		params.mBatchSize = B;
		params.mIterations = IT;
		params.mLearningRate = 1.0 / (1 << 3);

		const Decimal D = D16;
		aby3ML p;

		p.init(pIdx, chlPrev, chlNext, toBlock(pIdx));

		//p.mPrint = cmd.isSet("print");

		sf64Matrix<D> train_data, train_label, W2, test_data, test_label;

		if (pIdx == 0)
		{
			train_data = p.localInput<D>(val_train_data);
			train_label = p.localInput<D>(val_train_label);
			W2 = p.localInput<D>(val_W2);
			test_data = p.localInput<D>(val_test_data);
			test_label = p.localInput<D>(val_test_label);
		}
		else
		{
			train_data = p.remoteInput<D>(0);
			train_label = p.remoteInput<D>(0);
			W2 = p.remoteInput<D>(0);
			test_data = p.remoteInput<D>(0);
			test_label = p.remoteInput<D>(0);
		}


		p.mPreproNext.resetStats();
		p.mPreproPrev.resetStats();

		auto preStart = std::chrono::system_clock::now();

		p.preprocess((B + dim) * IT, D);

		double preBytes = p.mPreproNext.getTotalDataSent() + p.mPreproPrev.getTotalDataSent();


		p.mNext.resetStats();
		p.mPrev.resetStats();

		// printRevealedValues(W2, gen, print, dim, p);
		auto w2Val = p.reveal(W2);
		if (print) {
				for (u64 i = 0; i < dim; ++i) {
						std::cout << i << " " << gen.mModel(i, 0) << " " << w2Val(i, 0) << std::endl;
				}
		}

		auto start = std::chrono::system_clock::now();

		if (cmd.isSet("noOnline") == false)
			SGD_Logistic(params, p, train_data, train_label, W2, &test_data, &test_label, thre);
		//val_W2 = p.reveal(W2);

		//auto end = std::chrono::system_clock::now();


		//engine.sync();
		auto now = std::chrono::system_clock::now();
		auto preSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(start - preStart).count() / 1000.0;
		auto seconds = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count() / 1000.0;

		double bytes = p.mNext.getTotalDataSent() + p.mPrev.getTotalDataSent();

		if (print)
		{
			ostreamLock ooo(std::cout);
			ooo << "N: " << N << " D:" << dim << " B:" << B << " IT:" << IT << " => "
				<< (double(IT) / seconds) << "  iters/s  " << (bytes * 8 / 1024 / 2024) / seconds << " Mbps"
				<< " offline: " << (double(IT) / preSeconds) << "  iters/s  " << (preBytes * 8 / 1024 / 2024) / preSeconds << " Mbps" << std::endl;
		}


		auto w2Val = p.reveal(W2);

		if (print)
		{

			for (u64 i = 0; i < (u64)dim; ++i)
			{
				std::cout << i << " " << gen.mModel(i, 0) << " " << w2Val(i, 0) << std::endl;
			}
		}


		return 0;
	}


	int logistic_main_3pc_sh(oc::CLP & cmd)
	{

		auto N = cmd.getManyOr<int>("N", { 10000 });
		auto D = cmd.getManyOr<int>("D", { 1000 });
		auto B = cmd.getManyOr<int>("B", { 128 });
		auto IT = cmd.getManyOr<int>("I", { 10000 });
		auto testN = cmd.getOr("testN", 1000);

		IOService ios(cmd.isSet("p") ? 2 : 6);

		std::vector<std::thread> thrds;
		for (u64 i = 0; i < 3; ++i)
		{
			if (cmd.isSet("p") == false || cmd.get<u64>("p") == i)
			{
				thrds.emplace_back(std::thread([i, N, D, B, IT, testN, &cmd, &ios]() {

					auto next = (i + 1) % 3;
					auto prev = (i + 2) % 3;
					auto cNameNext = std::to_string(std::min(i, next)) + std::to_string(std::max(i, next));
					auto cNamePrev = std::to_string(std::min(i, prev)) + std::to_string(std::max(i, prev));

					auto modeNext = i < next ? SessionMode::Server : SessionMode::Client;
					auto modePrev = i < prev ? SessionMode::Server : SessionMode::Client;


					auto portNext = 1212 + std::min(i, next);
					auto portPrev = 1212 + std::min(i, prev);

					Session epNext(ios, "127.0.0.1", portNext, modeNext, cNameNext);
					Session epPrev(ios, "127.0.0.1", portPrev, modePrev, cNamePrev);

					std::cout << "party " << i << " next " << portNext << " mode=server?:" << (modeNext == SessionMode::Server) << " name " << cNameNext << std::endl;
					std::cout << "party " << i << " prev " << portPrev << " mode=server?:" << (modePrev == SessionMode::Server) << " name " << cNamePrev << std::endl;
					auto chlNext = epNext.addChannel();
					auto chlPrev = epPrev.addChannel();

					chlNext.waitForConnection();
					chlPrev.waitForConnection();

					chlNext.send(i);
					chlPrev.send(i);
					u64 prevAct, nextAct;
					chlNext.recv(nextAct);
					chlPrev.recv(prevAct);

					if (next != nextAct)
						std::cout << " bad next party idx, act: " << nextAct << " exp: " << next << std::endl;
					if (prev != prevAct)
						std::cout << " bad prev party idx, act: " << prevAct << " exp: " << prev << std::endl;

					ostreamLock(std::cout) << "party " << i << " start" << std::endl;

					auto print = cmd.isSet("p") || i == 0;

					for (auto n : N)
					{
						for (auto d : D)
						{
							for (auto b : B)
							{
								for (auto it : IT)
								{
									logistic_main_3pc_sh(n, d, b, it, testN, i, print, cmd, epPrev, epNext);
								}
							}
						}
					}
					}));
			}
		}

		for (auto& t : thrds)
			t.join();

		return 0;
	}

}

