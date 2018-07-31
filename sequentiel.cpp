#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#define dz 10
#define dt 1
#define T_initial_Bloc 825
#define T_initial_Croute 100
#define temps_Max 4500000// ~ 1 mois en second 2592000
#define Z 80000 //80000 m

const int N = Z / dz;
const int M = temps_Max / dt;

float  k = 14.19f, k1 = 28.38f, k2 = 26.27f;
class diffusionthermique
{
private:
	float** grille;
	int temps_actuel;
	bool is_stable();
	float calculNew_T(int);
public:
	diffusionthermique();
	void diffusion();	
};

diffusionthermique::diffusionthermique()
{
	grille = new float*[2];
	grille[0] = new float[N + 1];
	grille[1] = new float[N + 1];
	//initialisation de la grille et des temperature initial
	for (int m = 0; m<2; m++)
	{
		for (int z = 0; z <= N; z++)
		{
			grille[m][z] = 0.0f; 
		}
	}
	int z = 0;
	for (z; z*dz <= 11000; z++)
	{
		//initialisation temperature du bloc magnetique a  l'instant t=0
		grille[0][z] = 825.0f;
	}
	for (z; z*dz <= Z; z++)
	{
		//initialisation temperature de la surface et de la croute terestre a  l'instant t=0
		grille[0][z] = 100.0f;
	}
	grille[0][z] = 100.0f;
	temps_actuel = 0;
}

float diffusionthermique::calculNew_T(int z)
{
	float new_T;
	float K;

	if (z*dz <= 11000) //On est dans le bloc
	{
		if (grille[0][z] > 725)
		{
			K = k;
		}
		else
		{
			K = k1;
		}
	}
	else // On est dans la croute terteste
	{
		K = k2;
	}

	new_T = (grille[0][z] + K * ((grille[0][z + 1] - 2 * (grille[0][z]) + grille[0][z - 1]) / (dz*dz)));

	return new_T;
}

void diffusionthermique::diffusion()
{
	ofstream myfile("test.txt", ios::out);
	if (myfile.is_open())
	{
		bool first_iteration = true;
		do
		{
			if (first_iteration)
			{
				first_iteration = false;
			}
			else
			{
				for (int z = 1; z < N; z++) // Sauvegarde les nouvelles valeurs dans grille[0] pour calcul suivant
				{
					grille[0][z] = grille[1][z];
				}
			}


			temps_actuel++;

            
            if ((temps_actuel % 100000) == 0) // On ecrit dans le fichier a intervalle de 1000 iterations
            {
                cout << temps_actuel << endl;
                for (int z = 100; z*dz <= 80000 ; z+=100)
                {
                    myfile << setw(7) << grille[0][z] << " ";
                    
                    if ((z*dz)== 11000)
                    {
                        myfile << setw(7) << " ### "; // Met une barriere entre le block et la croute
                    }	
                }
                myfile << "\n\n";
            }
            
			for (int z = 1; z < N ; z++)
			{
				grille[1][z] = calculNew_T(z);
			}

		}while (!is_stable() && temps_actuel <= M);


		//AFFICHE L'ETAT DU BLOCK ET DE LA CROUTE STABILISE
		myfile << "L'ensemble s'est stabilise apres " << temps_actuel << " iterations" << endl << endl;
        
		for (int z = 100; z*dz <= 80000; z += 100)
		{
			myfile << setw(7) << grille[0][z] << " ";

			if ((z*dz)==11000)
			{
				myfile << setw(7) << " ### "; // Met une barriere entre le block et la croute
			}
		}
		myfile << "\n\n";

		cout << "L'ensemble s'est stabilise apres " << temps_actuel << " iterations" << endl << endl;
		myfile.close();
	}
	else
	{
		cout << "le Fichier n'a pas pu s'ouvrir..." << endl;
	}
}

//Verifie la stabilité au niveau du bloc
bool diffusionthermique::is_stable()
{
	bool stable = true;
	for (int position = 1; position*dz <= 11000; position++)
	{
		if (grille[0][position] != grille[1][position])
		{
			stable = false;
			break;
		}
	}
	return stable;
}

int main(int argc, char *argv[])
{
	diffusionthermique D;
	D.diffusion();
	return 0;
}