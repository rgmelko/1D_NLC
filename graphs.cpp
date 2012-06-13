#include "graphs.h"
#include <sstream>

graph::graph()
{
    NumberSites = -99;
    NumberBonds = -99;
    LatticeConstant = -99;
    Identifier = -99;
    RealSpaceCoordinates.clear();
    SubgraphList.clear();
    AdjacencyList.clear();
}

graph::graph(vector< pair<int, int> > & AdjList, int IdentNumber, int order, int edgeCount, 
	     int LattConst, vector< pair<int, int> > & subgraphs )
{
    AdjacencyList = AdjList;
    Identifier = IdentNumber;
    NumberSites = order;
    NumberBonds = edgeCount;
    LatticeConstant = LattConst;
    SubgraphList = subgraphs;
    RealSpaceCoordinates.clear();
}

graph::graph(vector< pair<int, int> > & AdjList, int IdentNumber, int order, int edgeCount, 
	     int LattConst, vector< pair<int, int> > & subgraphs, vector< vector< pair<int,int> > > embeddings )
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
    cout<<NumberSites<<" ";
    cout<<NumberBonds<<" ";
    cout<<LatticeConstant<<"\n";
    for (int i = 0; i< AdjacencyList.size(); i++){
        cout<<AdjacencyList[i].first<<" ";
        cout<<AdjacencyList[i].second<<"\n";
    }
    for (int i=0; i< SubgraphList.size(); i++){
      cout << SubgraphList[i].first<<" ";
      cout << SubgraphList[i].second<<"\n";
    }
      
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
            output<<"("<<graphList[currentGraph].AdjacencyList[currentBond].first<<","
                  <<graphList[currentGraph].AdjacencyList[currentBond].second<<")";
        }
        output<<endl;

        for (unsigned int currentSubgraph = 0; currentSubgraph < graphList[currentGraph].SubgraphList.size(); currentSubgraph++)
        {
	  //    output<<graphList[currentGraph].SubgraphList[currentSubgraph]<<" ";
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
    graph tempGraph;

    while ( !input.eof() )
    {
        rawLines.resize(rawLines.size() + 1);
        getline(input, rawLines.back()) ; 
    }

    input.close();

    //cout<<rawLines.size()<<endl;

    stringstream ss (stringstream::in | stringstream::out);

    for  (unsigned int currentLine = 0; currentLine < rawLines.size()-1; currentLine+=3)
      //for  (unsigned int currentLine = 0; currentLine < 9; currentLine+=3)
    {
        //cout<<currentLine<<" ";
        currentGraph = currentLine/memberCount;
        //cout<<currentGraph<<" ";
        unsigned int currentChar = 0;
        string currentNumber;
        //cout<<currentLine % memberCount<<endl;
    
        ss << rawLines.at(currentLine);
          
        ss >> tempGraph.Identifier;
        ss >> tempGraph.NumberSites;
        ss >> tempGraph.NumberBonds;
        ss >> tempGraph.LatticeConstant;
           
        //cout << "Identifier = " <<tempGraph.Identifier << endl;
        //cout << "NumberSites = " << tempGraph.NumberSites << endl;
        //cout << "Identifier = " <<tempGraph.NumberBonds << endl;
        //cout << "NumberSites = " << tempGraph.LatticeConstant << endl;

        ss.str("");
        ss.clear();

        //read in bonds

        
        ss << rawLines.at(currentLine+1);
        string teststring;
        tempGraph.AdjacencyList.resize(tempGraph.NumberBonds);
        for(int b=0; b<tempGraph.NumberBonds;b++){
          ss >> tempGraph.AdjacencyList[b].first;
          ss >> tempGraph.AdjacencyList[b].second;
          //cout << tempGraph.AdjacencyList[b].first << "," <<tempGraph.AdjacencyList[b].second << endl;
        }

        ss.str("");
        ss.clear();

        //read in subclusters
        ss << rawLines.at(currentLine+2);   
	int subSize(0);
	while(!ss.eof()){
	  ss >> teststring;
	  subSize++;
	}

	ss << rawLines.at(currentLine+2);
        tempGraph.SubgraphList.resize(5);
        for(int b=0; b<5;b++){
          ss >> tempGraph.SubgraphList[b].first;
          ss >> tempGraph.SubgraphList[b].second;
          //cout << tempGraph.AdjacencyList[b].first << "," <<tempGraph.AdjacencyList[b].second << endl;
        }

        ss.str("");
        ss.clear();

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


