/*
DATE : 15 AVRIL 2016
PROJET : IMPLEMENTION PARALLELE DU PROBLEME DE DIFFUSION DE CHALEUR AVEC MPI/C++
Evance Kafando
UNIVERSITE DE MONCTON
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <mpi.h>

using namespace std;

#define dz 10
#define dt 1
#define T_initial_Bloc 825
#define T_initial_Croute 100
#define temps_Max 4500000
#define Z 80000 

const int N = Z / dz;
const int M = temps_Max / dt;

//
int p, i, myid, reste , *nbrEnvoyer, *nbrToRecev, *deplacement_Envoi, *deplacement_Reception, iteration, root;
double wtime;

float  k = 14.19f, k1 = 28.38f, k2 = 26.27f;

float  grille [N + 1], /*Contient les donnees du bloc et de la croute */
	   grille_local0[N+1], /*Variable local a chaque processeur. les donnees locales pour le calcul*/          
	   grille_local1[N+1]; /*Variable local a chaque processeur. contient le resultat des calculs*/

/*Calcul les nouvelles temperatures. Prend en parametre une position locale, et le deplacement dans la grille par le processus correspondant permettant ainsi de verifier la position dans la grille*/
float calculNew_T(int, int); 


int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	
	ofstream myfile("test.txt", ios::out); /*Fichier d'ecriture des donnees*/

    root=0;
	iteration = 0;

	reste = (N+1)%p; /*Verifie si le nombre d'element est divisible par le nombre de processeur*/

    if (myid==root)
    {
        wtime = MPI::Wtime ( );
    }


    nbrEnvoyer = new int[p]; /*Tableau contenant le nombre d'element a envoyer a chaque processeur*/
    deplacement_Envoi = new int[p];/*Tableau contenant pour chaque processeur le deplacement a effectuer avant l'envoie des donnees*/
    nbrToRecev = new int[p]; /*Tableau contenant le nombre d'element que chaque processeur envoie au processeur root apres le calcul*/
    deplacement_Reception = new int[p]; /*Tableau contenant pour chaque processeur le deplacement a la reception des donnees par le processus root*/
  
    int deplac_E = 0; 
    int deplac_R = 1;
    /*Determine le nombre d'elements a envoyer a chaque processeur 
    et le nombre d'element que chaque processeur envoi au processeur apres le calcul*/
    for (i = 0; i < p; i++)
    {
        nbrEnvoyer[i] = ((N+1)/p)+1;

        if (i < reste) 
        {
            nbrEnvoyer[i] = nbrEnvoyer[i] + 1;
        }

        if ((reste != 0) && ((i != 0) && (i != (p-1))))
        {
        	 nbrEnvoyer[i] = nbrEnvoyer[i] + 1;
        }
        
        deplacement_Envoi[i] = deplac_E;       
        deplac_E = deplac_E + nbrEnvoyer[i]-2;

        nbrToRecev[i] = nbrEnvoyer[i]-2;
        deplacement_Reception[i] = deplac_R;
        deplac_R = deplac_R + nbrToRecev[i];
    }

    /*Root initialise les variables*/
    if (myid == root)
    {
		int z = 0;
		for (z; z*dz <= 11000; z++)
		{
			grille[z] = 825.0f;
		}

		for (z; z*dz <= Z; z++)
		{
			grille[z] = 100.0f;
		}

		grille[z] = 100.0f;
		grille[z+1] = 100.0f;
    }

    /*Tant que le temps maximum n'est pas atteint ou que le block n'est pas stable*/
    while (iteration <= temps_Max)
    {
    	iteration++;
    	/*Envoie collective des donnes a chaque processeur selon la taille des donnees correspondant0*/
		MPI_Scatterv(&grille, nbrEnvoyer, deplacement_Envoi, MPI_FLOAT, &grille_local0, nbrEnvoyer[myid],MPI_FLOAT,root, MPI_COMM_WORLD);

		/*Chaque processeur effectue les calcul sur ses donnees*/
		for (int z = 1; z < (nbrEnvoyer[myid]-1) ; z++) 
			{
				grille_local1[z] = calculNew_T(z,deplacement_Envoi[myid]);
			}

		/*Chaque processeur envoie ces resultat au processus root*/	
		MPI_Gatherv(&grille_local1[1],nbrToRecev[myid], MPI_FLOAT, &grille, nbrToRecev, deplacement_Reception, MPI_FLOAT, root, MPI_COMM_WORLD); 

		/*Le Processus root ecrit les resulats a chaque 100000 ième iteration*/
		if ((myid == root) && ((iteration % 100000) == 0)) 
			{
				cout << iteration << endl;
				myfile << "0 ";
				for (int z = 100; z*dz <= 80000 ; z+=100) 
				{
					myfile << setw(7) << grille[z] << " ";

					if ((z*dz)== 11000)
					{
						myfile << setw(7) << " ### "; 
					}	
				}
				myfile << setw(7) << grille[N] << " ";
				myfile << "\n\n";
			}	
    }

    /*Le processus root affiche et ecrit les resultats dans le fichier*/
    if (myid == root)
    {
		myfile << "L'ensemble s'est stabilise apres " << iteration << " iterations" << endl << endl;
	
		myfile << grille[0] << " ";
		for (int z = 1000; z*dz <= 80000; z += 1000)
		{
			myfile << setw(7) << grille[z] << " ";

			if (((z*dz) % 11000) == 0)
			{
				myfile << setw(7) << " ### ";
			}
		}
		myfile << setw(7) << grille[N] << " ";
		myfile << "\n\n";

    	cout << "L'ensemble s'est stabilise apres " << iteration << " iterations" << endl << endl;
		myfile.close();

        wtime = MPI::Wtime ( ) - wtime;
        
        cout << "\n";
        cout << "  Wall clock elapsed seconds = " << wtime << "\n";
    }

	MPI_Finalize();
	return 0;
}

/*Calcul les nouvelles temperatures. Prend en parametre une position locale, et le deplacement dans la grille par le processus correspondant permettant ainsi de verifier la position dans la grille*/
float calculNew_T(int z, int deplacement_Envoi)
{
	float new_T;
	float K;
	/*Determine le facteur de diffusion approprié*/
	if ((z+deplacement_Envoi)*dz <= 11000)
	{
		if (grille_local0[z] > 725)
		{
			K = k;
		}
		else
		{
			K = k1;
		}
	}
	else 
	{
		K = k2;
	}
	/*Formule de la diffusion*/
	new_T = (grille_local0[z] + K * ((grille_local0[z + 1] - 2 * (grille_local0[z]) + grille_local0[z - 1]) / (dz*dz)));

	return new_T;
}

/*Verifie la stabilité au niveau du bloc
bool is_stable()
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
}*/