/*
 * dummy.cpp
 *
 *  Created on: Jul 22, 2022
 *      Author: brandon
 */


#include <commands/command.h>

struct CommandDummy : public Command
{
	CommandDummy() : Command("dummy", "jdftx/Electronic/Optimization")
	{
		format = "";
		comments = "";
	}

	void process(ParamList& pl, Everything& e)
	{
		std::cout <<"Hi!!";
	}

	void printStatus(Everything& e, int iRep)
	{
	}
}
commandDummy;
