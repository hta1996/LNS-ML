#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include "LNS.h"
#include "AnytimeBCBS.h"
#include "AnytimeEECBS.h"
#include "PIBT/pibt.h"
#include <fstream>
#include <iostream> 
using namespace std;    
//#include <torch/script.h>

/* Main function */
int main(int argc, char** argv)
{
    
    //torch::jit::script::Module module;
    ///module = torch::jit::load("traced_resnet_model.pt");
    //cerr<<"start"<<endl;
	namespace po = boost::program_options; 
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
        ("model",po::value<string>()->default_value("NONE"),"path to model file")
		// params for the input instance and experiment settings
		("map,m", po::value<string>()->required(), "input file for map")
		("agents,a", po::value<string>()->required(), "input file for agents")
		("agentNum,k", po::value<int>()->default_value(0), "number of agents")
        ("output,o", po::value<string>(), "output file")
		("cutoffTime,t", po::value<double>()->default_value(7200), "cutoff time (seconds)")
		("screen,s", po::value<int>()->default_value(0),
		        "screen option (0: none; 1: LNS results; 2:LNS detailed results; 3: MAPF detailed results)")
		("stats", po::value<string>(), "output stats file")

		// solver
		("solver", po::value<string>()->default_value("LNS"), "solver (LNS, A-BCBS, A-EECBS)")

        // params for LNS
        ("ubns", po::value<int>()->default_value(5), "upper bound of the neighborhood size")
        ("lbns", po::value<int>()->default_value(5), "lower bound of the neighborhood size")

        ("neighborSize", po::value<int>()->default_value(5), "Size of the neighborhood")
        ("maxIterations", po::value<int>()->default_value(1000000), "maximum number of iterations")
        ("initAlgo", po::value<string>()->default_value("EECBS"),
                "MAPF algorithm for finding the initial solution (EECBS, CBS, PP)")
        ("replanAlgo", po::value<string>()->default_value("CBS"),
                "MAPF algorithm for replanning (EECBS, CBS, PP)")
        ("destoryStrategy", po::value<string>()->default_value("RandomWalk"),
                "Heuristics for finding subgroups (RandomWalk, Intersection, Adaptive)")
        ("pibtWindow", po::value<int>()->default_value(5),
             "window size for winPIBT")
        ("winPibtSoftmode", po::value<bool>()->default_value(true),
             "winPIBT soft mode")
        ("numModels", po::value<int>()->default_value(1), "#model in the portfolio")
        ("featureFile",po::value<string>()->default_value("NONE"), "featureFile")
        ("rand",po::value<int>()->default_value(0),"random seed for agent permutation (must be a prime number larger than 1000), or 0 for non-randomization")
        ("initSol",po::value<string>()->default_value("NONE"),"file of initial solution")
        ("outputSol",po::value<string>()->default_value("NONE"),"file to output solution")
        ("samples", po::value<int>()->default_value(1), "#samples when collecting data")
        ("testing", po::value<int>()->default_value(0), "testing or not")
        ("port", po::value<int>()->default_value(0), "port number")
        ("weight", po::value<string>()->default_value("NONE"), "file of the model weight")

		;
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if (vm.count("help")) {
		cout << desc << endl;
		return 1;
	}

    PIBTPPS_option pipp_option;
    pipp_option.windowSize = vm["pibtWindow"].as<int>();
    pipp_option.winPIBTSoft = vm["winPibtSoftmode"].as<bool>();

    po::notify(vm);

	if(vm["rand"].as<int>()==0)srand((int)time(0));
    else srand(vm["rand"].as<int>());
    
    //cerr<<"loading instance"<<endl;

	Instance instance(vm["map"].as<string>(), vm["agents"].as<string>(), vm["rand"].as<int>(),
		vm["agentNum"].as<int>());
    double time_limit = vm["cutoffTime"].as<double>();
    int screen = vm["screen"].as<int>();
	//srand(0);
    //cerr<<"loaded instance"<<endl;
	if (vm["solver"].as<string>() == "LNS")
    {
        LNS lns(instance, time_limit,
                vm["initAlgo"].as<string>(),
                vm["replanAlgo"].as<string>(),
                vm["destoryStrategy"].as<string>(),
                vm["neighborSize"].as<int>(),
                vm["maxIterations"].as<int>(), screen, pipp_option);
        lns.num_of_samples=vm["samples"].as<int>();
        lns.num_of_model=vm["numModels"].as<int>();
        lns.testing=vm["testing"].as<int>();
        lns.PORT=vm["port"].as<int>();
        //cerr<<vm["weight"].as<string>()<<endl;
        if (vm["weight"].as<string>()!="NONE")
        {
            memset(lns.model_weight,0,sizeof(lns.model_weight));
           
            ifstream myfile(vm["weight"].as<string>().c_str());
			while(!myfile.eof())
			{
                int id; double wei;
				if(scanf("%d%lf",&id,&wei)==EOF)break;
				lns.model_weight[id-1]=wei;
                cerr<<id<<" "<<wei<<" ";
			}
			myfile.close();

           /*freopen(vm["weight"].as<string>().c_str(),"r",stdin);
            int id; double wei;
            while(scanf("%d%lf",&id,&wei)!=EOF)
            {
                lns.model_weight[id-1]=wei;
                cerr<<id<<" "<<wei<<" ";
            }
            fclose(stdin);*/
        }
        //if(lns.num_of_samples==1)
        {
            lns.ub_neighborsize=vm["neighborSize"].as<int>();
            lns.lb_neighborsize=vm["neighborSize"].as<int>();
        }//else
        if(vm["ubns"].as<int>()>vm["lbns"].as<int>())
        {
            lns.ub_neighborsize=vm["ubns"].as<int>();
            lns.lb_neighborsize=vm["lbns"].as<int>();
        }
        //cerr<<lns.num_of_samples<<endl;
        lns.model=vm["model"].as<string>();
        lns.feature_file=vm["featureFile"].as<string>();
        lns.initSol_file=vm["initSol"].as<string>();
        lns.outputSol_file=vm["outputSol"].as<string>();
        freopen(lns.feature_file.c_str(),"w",stdout);
        fclose(stdout);
        bool succ = lns.run(); 
        //if (succ)
          //  lns.validateSolution();
        if (vm.count("output"))
            lns.writeResultToFile(vm["output"].as<string>());
        if (vm.count("stats"))
            lns.writeIterStatsToFile(vm["stats"].as<string>());
    }
    else if (vm["solver"].as<string>() == "A-BCBS") // anytime BCBS(w, 1)
    {
        AnytimeBCBS bcbs(instance, time_limit, screen);
        bcbs.run();
        bcbs.validateSolution();
        if (vm.count("output"))
            bcbs.writeResultToFile(vm["output"].as<string>());
        if (vm.count("stats"))
            bcbs.writeIterStatsToFile(vm["stats"].as<string>());
    }
    else if (vm["solver"].as<string>() == "A-EECBS") // anytime EECBS
    {
        AnytimeEECBS eecbs(instance, time_limit, screen);
        eecbs.run();
        eecbs.validateSolution();
        if (vm.count("output"))
            eecbs.writeResultToFile(vm["output"].as<string>());
        if (vm.count("stats"))
            eecbs.writeIterStatsToFile(vm["stats"].as<string>());
    }
	else
    {
	    cerr << "Solver " << vm["solver"].as<string>() << " does not exist!" << endl;
	    exit(-1);
    }
	return 0;

}