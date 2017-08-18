#pragma once
#define SimCase_h

#include "BeamCSD.h"




class Cant :public ModelCase
{
public:
	Cant() :ModelCase() { ; }
	~Cant() { ; }
	Cant(const Cant &C) :ModelCase(C) { ; }

	void GetModel(void);
	void GetProb(void);

};

