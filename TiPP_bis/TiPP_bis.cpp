#include <iostream>
#include <fstream>
#include <string>
#include "class.cpp"
#define M_PI  3.14159265358979323846  /* pi */
#include <fstream>

// <math.h>
int main()
{

    //Init 
    std::cout << "Select the nuclei Sc:0 Ne: 1 other 2" << std::endl;
    int choix = 0; std::cin >> choix;
    std::cout << "Param  predefinis ? ->0" << std::endl;
    int choix2 = 0; std::cin >> choix2;
    double deltaE = 0.2; int Jmax = 40; double Emax = 40.0; double f_const = 0; std::string radical; std::string file_add_level_pos; std::string file_add_level_neg; std::string file_add_trans;
    if (choix == 2) 
    {
        std::cout << "enter the name of the level pos file "; std::cin >> file_add_level_pos ;
        std::cout << "enter the name of the level neg file "; std::cin >> file_add_level_neg;

        std::cout << "enter the name of the level transition "; std::cin >> file_add_trans; radical = file_add_trans;

    }
    if (choix2 !=0) 
    {
        std::cout << "Emax ?\n"; std::cin >> Emax; std::cout << "deltaE ?\n"; std::cin >> deltaE; std::cout << "Jmax ?\n"; std::cin >> Jmax;
        if (choix == 2) { std::cout << "f_const ?\n"; std::cin >> f_const; }
    }
    if (choix == 0) radical = "Sc44";
    else  radical = "Ne26";

    //creating of the first class
    rho R(Emax, deltaE, Jmax);



    // Start extracting the info from the documents
    double Emin= 0;  
    double temp; double temp2; int J; std::string line;
    std::ifstream file2;
    if(choix==0) file2.open("input/SC44POS");   //("C:\Users\seign\Desktop\Admin\TIPP\SC44POS"); 
    else if (choix == 1) file2.open("input/NE26POS_SPEC");
    else  file2.open("input/"+file_add_level_pos);
    while (!(file2.eof()) )                                                                                                                          ////////Read first doc
    {
        file2.ignore(256, '=');file2.ignore(256, '=');

        file2 >> J; //std::cout << J << " ";
        file2.ignore(256, '='); 
        file2.ignore(256, '='); file2.ignore(256, '='); file2.ignore(256, '=');
        file2 >> temp;
        file2.ignore(256, '=');
        file2 >> temp2; if (Emin >temp2)Emin = temp2;
        std::getline(file2, line);
        R.Add(J, temp);

    }
    file2.close();
    double Emin2=0;
    std::ifstream file;
    if(choix!=1)
    {
        if (choix == 0)file.open("input/SC44NEG");//(\SC44NEG");
        else  file.open("input/" + file_add_level_neg);
        while (!(file.eof()))
        {
            file.ignore(256, '='); file.ignore(256, '=');
            file >> J; //std::cout << J << std::endl;
            file.ignore(256, '='); file.ignore(256, '='); file.ignore(256, '='); file.ignore(256, '=');
            file >> temp;
            file2.ignore(256, '=');
            file2 >> temp2;
            if (Emin2 > temp2)Emin2 = temp2;
            std::getline(file, line);
            R.Add(J, temp - (Emin2 - Emin));

        }

    }
    file.close();
    //Display of rho
    int size;
    std::cout << "total number of level recorded: " << R.total() << std::endl;
    std::ofstream myfile;
    myfile.open("output/rho_" + radical + ".txt.txt");

    size = (int)(Emax / deltaE) + 1;
    for (int i = 0; i < size; i++)
    {
        //std::cout << R.SumJ(i) << " ";
        myfile << i * deltaE << "; " << R.SumJ(i) << std::endl;
    }
    myfile.close();



    // reading the transition file
    {
        int J1; int J2; double E; double dE; double redprob1; double redprob2;
        std::ifstream file3;
        if(choix==0)  file3.open("input/sc44e1");//("\sc44e1");   
        else if (choix == 1) file3.open("input/ne26pos");
        else  file3.open("input/"+file_add_trans);

        B bb(Emax, deltaE, Jmax,Emin);
        int error = 0;
        while (!(file3.eof()))
        {
            file3 >> J1;
            file3 >> J2;
            file3 >> E;
            file3 >> dE;
            file3 >> redprob1;
            file3 >> redprob2;
            std::getline(file3, line);
            //std::cout << J1 << " " << J2 << " " << dE << std::endl;

            if (dE > 0)
            {
                error += bb.Add(J1, dE, redprob1, E - Emin); // error sum the time the operation add was not possible

                //std::cout << dE << " ";
            }
            else
                error += bb.Add(J2, -1 * dE, redprob2, E - dE - Emin);
        }
        std::cout <<"number of transition outside of the predetermied parameters: "<< error << std::endl;
        //bb.Avg();

        std::ofstream myfile2;
        myfile2.open("output/B_"+ radical+".txt");
        std::ofstream myfile3;
        myfile3.open("output/f_" + radical + ".txt");
        double f_const = 0; int c = 0; double res; int total = 0;
          if(choix==0)  f_const = 137  * M_PI / (pow(9 * 197.3269805, 3)) * 16;
          if (choix == 1)f_const = 1.154 * pow(10,-8);
        for (int i = 0; i < size; i++)
        {
            //std::cout << bb.read(1, i) << std::endl;
            myfile2 << i * deltaE << "; " << bb.Avg(i) /*<< "; " << bb.Avg(i,0.0,5.0) << "; " << bb.Avg(i,5.0,10) << "; " << bb.Avg(i,10.0,15) */<< std::endl;      ////////Write <B(E1)>
            //std::cout << 137 * 16 / (pow(9 * 197.3269805, 3)) << std::endl;
            res = f_const * bb.f(i, R, c);
            total += c;
            myfile3 << i * deltaE << "; " <<res <<" ; "<<c<< std::endl;
       }
        std::cout <<"f const : "<< f_const <<" ; total number of transition considered: "<<total<< std::endl;
        myfile2.close();
    }
    std::cout << "EXIT" << std::endl;
    std::cin >> Emax; //command to finsih the programme
}