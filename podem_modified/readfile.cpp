	#include <iostream>
#include <string>
#include <vector>
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
#include <sstream>
#include <vector>
#include <random>
using namespace std;

int main(int argc, char const *argv[])
{
vector<string> P2_words;
	ifstream Phase2("c1196_opt phase II");

	vector<vector<char> > P2_Final_Indiv;

    string word;
	while (Phase2 >> word)
	{
		P2_words.push_back(word);
	}
	for(int s = 14003; s < P2_words.size(); s++){
			vector<char> Final_gene(P2_words[s].begin(), P2_words[s].end());
			if (s == 14003)
			{
				Final_gene.erase(Final_gene.begin());
			}
			P2_Final_Indiv.push_back(Final_gene);
	}
	cout << P2_Final_Indiv.size()<<endl;
    for(int m = 0; m < P2_Final_Indiv.size(); m++){
        for(int n = 0; n < P2_Final_Indiv[m].size(); n++){
            cout <<P2_Final_Indiv[m][n] ;
        }

        cout << endl;
    }
		/* code */
	return 0;
}
	