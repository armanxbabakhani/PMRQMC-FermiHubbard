#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>

using namespace std;

struct FHdata{
    double t , U;
    vector<vector<bool>> Adjacency;
};



void Print_offdiag(ofstream& output , int i , int j , double t){
    // IMPORTANT NOTE: here it is implicitly assumed that i < j
    int tail = pow(-1 , j - i + 1);
    string minusHalft = to_string(-0.5*t*tail) , plusHalft = to_string(0.5*tail);
    output << minusHalft << " ";
    for(int delta = 0; delta < j - i + 1; delta++){
        output << to_string(i + delta) << " 3 1 2 ";
    }
    output << i << " 1 1 2 " << j << " 1 1 2"<< endl;
    output << plusHalft << " ";
    for(int delta = 1; delta < j-i; delta++){
        output << to_string(i + delta) << " 3 1 2 ";
    }
    output << i << " 1 1 2 " << j << " 1 1 2" << endl;
}

FHdata Data_extract(const string& filename) {
    ifstream inputFile(filename);
    FHdata AllData;
    vector<pair<int, vector<int>>> AdjacencyList;

    if (!inputFile.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        return AllData;
    }

    // Read the first two lines to extract t and U
    double t , U;
    int NodeCount = 0;
    string line;
    while(getline(inputFile, line)){
        istringstream iss(line);
        vector<string> tokens(std::istream_iterator<string>{iss},std::istream_iterator<string> ());
        if(tokens[0][0] == 't'){
            t = stod(tokens[2]);
        }
        else if(tokens[0][0] == 'U'){
            U = stod(tokens[2]);
        }
        // Counting the nodes on the graph
        else{
            //size_t ColonPos = line.find(':');
            //string Node = line.substr(0, ColonPos);
            int CurrentNode , node , nodestart;
            bool coloncontained = true;
            if(tokens[0][tokens[0].size()-1] == ':'){
                CurrentNode = stoi(tokens[0].substr(0,tokens[0].size()-1));
            }
            else{
                CurrentNode = stoi(tokens[0]);
                coloncontained = false;
            }
            pair<int, vector<int>> NodeList;
            vector<int> adjacency;
            if(CurrentNode > NodeCount){
                NodeCount = CurrentNode;
            }
            NodeList.first = CurrentNode;

            if(coloncontained){
                nodestart = 1;
            }
            else{
                nodestart = 2;
            }
            for(int i=nodestart; i < tokens.size();i++){
                adjacency.push_back(stoi(tokens[i]));
            }
            NodeList.second = adjacency;
            AdjacencyList.push_back(NodeList);
        }
    }

    // Building the Adjacency matrix now
    vector<vector<bool>> AdjacencyMatrix(NodeCount , vector<bool>(NodeCount , 0));
    for(int i = 0; i < AdjacencyList.size(); i++){
        pair<int , vector<int>> CurrentNodeList = AdjacencyList[i];
        int CurrentNode = CurrentNodeList.first;
        vector<int> NodeList = CurrentNodeList.second;
        for(int j = 0; j < NodeList.size(); j++){
            AdjacencyMatrix[CurrentNode-1][NodeList[j]-1] = 1;
        }
    }
    inputFile.close();

    AllData.t = t;
    AllData.U = U;
    AllData.Adjacency = AdjacencyMatrix;
    return AllData;
}


void FH_PMR_convert(FHdata AllData , string filename){
    double t = AllData.t , U = AllData.U;
    vector<vector<bool>> A = AllData.Adjacency;
    int N = A.size();
    cout << "Number of particles is " << N << endl;
    string Uover4 = to_string(U/4) , UNover4 = to_string(U*N/4);

    ofstream OutputFile(filename);
    // Print the Diagonal part
    if(!OutputFile.is_open()){
        cerr << "Error opening the file: " << filename << endl;
    }

    // Havent included the constant shift of U*N/4 (yet)
    if(abs(U)>1E-6){
        OutputFile << UNover4 << " 1 0 1 2" << endl;  
        // The diagonal (D0) part:
        for(int i = 1; i < N+1; i++){
            OutputFile << Uover4 << " " << i << " 3 1 2" << endl;
            OutputFile << Uover4 << " " << i + N << " 3 1 2" << endl;
            OutputFile << Uover4 << " " << i << " 3 1 2 " << i + N << " 3 1 2" << endl; 
        }
    }
    
    // The off-diagonal part:
    for(int i = 0; i < N; i++){
        for(int j = i+1; j < N; j++){
            if(A[i][j]){
                Print_offdiag(OutputFile , i + 1, j + 1 , t); // For spin up
                Print_offdiag(OutputFile , i + N + 1 , j + N + 1 , t); // For spin down
            }
        }
    } 
    OutputFile.close();
}

int main(int argc , char* argv[]){
    string filename(argv[1]);  // Reading the name of the input .txt file describing the Hamiltonian
    FHdata Data;
    Data = Data_extract(filename);

    size_t pos = filename.rfind('.');
    string output = filename.substr(0, pos);

    string OutputName(output+"_PMRQMC.txt");
    FH_PMR_convert(Data , OutputName);
    return 0;
    
}