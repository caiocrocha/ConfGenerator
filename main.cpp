#include <iostream>
#include <fstream>
#include <iomanip>
#define MAX 200

using namespace std;

bool letra(string str)
{
    return ((str[0] >= 'A' && str[0] <= 'Z')||(str[0] >= 'a' && str[0] <= 'a'));
}

bool alg(string str)
{
    return ((str[0] >= '0' && str[0] <= '9')
            ||(str[0] == '-' && (str[1] >= '0' && str[1] <= '9')));
}

int main ()
{
    double vet[MAX];
    string atual, str;
    ifstream input("butane.txt");
    fstream output("saida.txt", fstream::app);
    if(input.is_open() && output.is_open())
    {
        int n = 0;
        for(int i = 0; i < MAX && input >> atual; i++)
        {
            if(letra(atual))
            {
                str = atual;
                output << str << endl;
            }
            else if(alg(atual))
            {
                vet[n] = atof(atual.c_str());
                output << fixed << setprecision(6) << vet[n] << endl;
                n++;
            }
        }
        input.close();
        output.close();
        cout << str << endl;
        for(int i = 0; i < n; i++)
            cout << vet[i] << " ";
        cout << endl;
    }
    else
        cout << "Unable to open files" << endl;
    /*
    int n = 10, p = 20;
    int x[n][p];
    ifstream myfile;
    ofstream output("saida.txt");
    myfile.open("numeros.txt");
    string line;
    string temp = "";
    int a = 0; // row index
    if(output.is_open()){
    while (getline(myfile, line))   //while there is a line
    {
        int b = 0; // column index
        for (unsigned int i = 0; i < line.size(); i++)   // for each character in rowstring
        {
            if (!isspace(line[i]))   // if it is not blank, do this
            {
                string d(1, line[i]); // convert character to string
                temp.append(d); // append the two strings
            }
            else
            {
                x[a][b] = stod(temp);  // convert string to double
                temp = ""; // reset the capture
                b++; // increment b cause we have a new number
                output << x[a][b];
            }
        }
        x[a][b] = stod(temp);
        temp = "";
        a++; // onto next row
    }}
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < p; j++)
            cout << x[i][j] << " ";
        cout << endl;
    }
    */
    return 0;
}
