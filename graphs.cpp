#include "graphs.h"

graph::graph()
{
	NumberSites = 0;
	NumberBonds = 0;
	LatticeConstant = 1;
	Identifier = 0;
	RealSpaceCoordinates.clear();
	SubgraphList.clear();
	AdjacencyList.clear();
}

graph::graph(vector< pair<int, int> > & AdjList, int IdentNumber, int order, int edgeCount, int LattConst, vector< int > & subgraphs )
{
	AdjacencyList = AdjList;
	Identifier = IdentNumber;
	NumberSites = order;
	NumberBonds = edgeCount;
	LatticeConstant = LattConst;
	SubgraphList = subgraphs;
	RealSpaceCoordinates.clear();
}

graph::graph(vector< pair<int, int> > & AdjList, int IdentNumber, int order, int edgeCount, int LattConst, vector< int > & subgraphs, vector< vector< pair<int,int> > > embeddings )
{
	AdjacencyList = AdjList;
	Identifier = IdentNumber;
	NumberSites = order;
	NumberBonds = edgeCount;
	LatticeConstant = LattConst;
	SubgraphList = subgraphs;
	RealSpaceCoordinates = embeddings;
}

void graph::print()
{
    cout<<Identifier<<" ";
    cout<<NumberBonds<<" ";
    cout<<NumberSites<<" ";
    cout<<LatticeConstant<<"\n";
};


graph& graph::operator=( const graph & other)
{
    this->NumberBonds = other.NumberBonds;
    this->NumberSites = other.NumberSites;
    this->LatticeConstant = other.LatticeConstant;
    this->AdjacencyList = other.AdjacencyList;
    this->SubgraphList = other.SubgraphList;
    return *this;
};

bool graph::operator==( const graph & other)
{
	return (( this->NumberBonds == other.NumberBonds) &&
			( this->NumberSites == other.NumberSites) &&
			( this->LatticeConstant == other.LatticeConstant) &&
			( this->AdjacencyList == other.AdjacencyList) );
}


void WriteGraphsToFile( vector< graph > & graphList, std::string file)
{
	ofstream output(file.c_str());
	for( unsigned int currentGraph = 0; currentGraph < graphList.size(); currentGraph++)
	{
		output<<graphList[currentGraph].Identifier<<endl;
		output<<graphList[currentGraph].NumberSites<<endl;
		output<<graphList[currentGraph].NumberBonds<<endl;
		output<<graphList[currentGraph].LatticeConstant<<endl;

		for (unsigned int currentBond = 0; currentBond < graphList[currentGraph].AdjacencyList.size(); currentBond++)
		{
			output<<"("<<graphList[currentGraph].AdjacencyList[currentBond].first<<","<<graphList[currentGraph].AdjacencyList[currentBond].second<<")";
		}
		output<<endl;

		for (unsigned int currentSubgraph = 0; currentSubgraph < graphList[currentGraph].SubgraphList.size(); currentSubgraph++)
		{
			output<<"("<<graphList[currentGraph].SubgraphList[currentSubgraph]<<")";
		}
		output<<endl;
	}
}

void ReadGraphsFromFile( vector< graph > & graphList, const string & file)
{
	ifstream input(file.c_str());
	vector< string > rawLines;
	int currentGraph;
	const int memberCount = 6;

	while ( !input.eof() )
	{
		rawLines.resize(rawLines.size() + 1);
		getline(input, rawLines.back()) ; 
	}

	input.close();

	//cout<<rawLines.size()<<endl;

	for (unsigned int currentLine = 0; currentLine < rawLines.size(); currentLine++)
	{
		cout<<currentLine<<" ";
		currentGraph = currentLine/memberCount;
		cout<<currentGraph<<" ";
		graph tempGraph;
		unsigned int currentChar = 0;
		string currentNumber;
	    cout<<currentLine % memberCount<<endl;	
		switch ( currentLine % memberCount )
		{
			case 0 :
				tempGraph.Identifier = atoi(rawLines.at(currentLine).c_str());
			case 1 :
				tempGraph.NumberSites = atoi(rawLines.at(currentLine).c_str());
				break;
			case 2 :
				tempGraph.NumberBonds = atoi(rawLines.at(currentLine).c_str());
				break;
			case 3 : 
				tempGraph.LatticeConstant = atoi(rawLines.at(currentLine).c_str());
				break;
			case 4 : 

				while ( currentChar < rawLines.at(currentLine).length() )
				{
					if ( rawLines.at(currentLine)[currentChar] == '(' )
					{
						tempGraph.AdjacencyList.resize( tempGraph.AdjacencyList.size() + 1);
					}
					if ( rawLines.at(currentLine)[currentChar] != '(' &&
							rawLines.at(currentLine)[currentChar] != ')' &&
							rawLines.at(currentLine)[currentChar] != ',' && 
							rawLines.at(currentLine)[currentChar] != '\n' )
					{
						currentNumber.push_back(rawLines.at(currentLine)[currentChar]);
					}
					if ( rawLines.at(currentLine)[currentChar] == ',' )
					{
						tempGraph.AdjacencyList.back().first = atoi(currentNumber.c_str());
						currentNumber.clear();
					}
					if ( rawLines.at(currentLine)[currentChar] == ')' )
					{
						tempGraph.AdjacencyList.back().second = atoi(currentNumber.c_str());
						currentNumber.clear();
					}
				}
				break;
			case 5 :

				while ( currentChar < rawLines.at(currentLine).length() )
				{
					if ( rawLines.at(currentLine)[currentChar] == '(' )
					{
						tempGraph.SubgraphList.resize( tempGraph.SubgraphList.size() + 1);
					}

					if ( rawLines.at(currentLine)[currentChar] != '(' &&
							rawLines.at(currentLine)[currentChar] != ')' &&
							rawLines.at(currentLine)[currentChar] != '\n' )
					{
						currentNumber.push_back(rawLines.at(currentLine)[currentChar]);
					}
					if ( rawLines.at(currentLine)[currentChar] == ')' )
					{
						tempGraph.SubgraphList.back() = atoi(currentNumber.c_str());
						currentNumber.clear();
					}
				}
				break;
		}

		graphList.push_back(tempGraph);

	}
}


graph GetGraphFromFile(const int IdNumber, const string & file)
{
	vector< graph > fileGraphs;
	ReadGraphsFromFile( fileGraphs, file);
	for( unsigned int currentGraph = 0; currentGraph < fileGraphs.size(); currentGraph++)
	{
		if ( fileGraphs.at(currentGraph).Identifier == IdNumber)
		{
			return fileGraphs.at(currentGraph);
		}
	}
}


