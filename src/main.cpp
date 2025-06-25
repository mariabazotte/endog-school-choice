#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <iterator>  

#include "input/input.hpp"
#include "instance/instance.hpp"
#include "algorithm/algofactory.hpp"
#include "algorithm/algorithm.hpp"

int main(int argc, char *argv[]) {
	/* Data information. */
	try{
		Input input(argc,argv);
		input.display();
        Instance instance(input);
		instance.display();
		AlgoFactory factory;
		Algorithm * algo = factory.createAlgorithm(input,instance);
        algo->solve();
		input.writeHead(algo->writeHead());
		input.write(algo->write());
		delete algo;
	}catch (GRBException e){
    	std::cout << "Error code = " << e.getErrorCode() << std::endl;
    	std::cout << e.getMessage() << std::endl;
  	}
	catch (const std::string & e) { std::cout << "EXCEPTION | " << e << std::endl; }
	catch (const std::exception & e) { std::cout << "EXCEPTION | " << e.what() << std::endl; }

    return 0; 
}