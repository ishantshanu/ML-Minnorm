#pragma once

class MyEntryPoint
{
public:
	static void MyMain(std::vector< std::vector< int > >cliqueinfo, int total_nodes, std::vector< double >x_neg_V);

private:

	void WriteNumOnesHeader();
	void StereoScriptTest();
	void StereoFusionMoveScriptTest();
	void FOEScriptTest();
	void DenoiseAlphaExpansionScript();
};

