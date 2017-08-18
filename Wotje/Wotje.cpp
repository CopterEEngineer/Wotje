// Wotje.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "BeamCSD.h"
#include "SimCase.h"

int main()
{
	Beam_1D beam1d;
	Cant cant;

	cant.GetModel();
	cant.GetProb();
	beam1d.InitBeam1D(cant);
	beam1d.AssembleKM();
	if (beam1d.SetBdryCond())
	{
		beam1d.StaticSolution();
		beam1d.ModeSolution();
		beam1d.DynamicSolution();
	}
	else
		printf("Boundary condition setting unsuccessfully. \n");

	
	system("pause");
	return 0;
}

