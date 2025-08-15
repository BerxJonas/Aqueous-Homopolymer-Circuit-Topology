#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <err.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>
//#include <omp.h>

void multichain_topology(double bl, int nr_contacts, int nr_beads, int nr_polymers, int nr_atoms, int contacts[((int) nr_atoms*(nr_atoms-1))/2][4], double coordinates[nr_beads][3][nr_polymers],int *ct, int *ls);
void contact_map(FILE *fp, int nr_contacts, int nr_beads, int nr_polymers, int nr_atoms, int contacts[((int) nr_atoms*(nr_atoms-1))/2][4], double coordinates[nr_beads][3][nr_polymers], int *ct, int *ls);
int motif(int poly1, int poly2, int poly3, int poly4, int pos1, int pos2, int pos3, int pos4);
int loop_size(int nr_beads, int nr_polymers, int poly1, int poly2, int poly3, int poly4, int pos1, int pos2, int pos3, int pos4, int ct_motif);
double bond_length(int nr_beads, int nr_polymers, double coordinates[nr_beads][3][nr_polymers]);

int main(int argc, char *argv[]){
	clock_t start, end;
     	double cpu_time_used;
     	start = clock();

	int nr_polymers = atoi(argv[1]); // # Polymers
	int nr_beads = atoi(argv[2]); // # Atoms per polymer
	int nr_atoms = nr_beads*nr_polymers;
	double cutoff = strtod(argv[3],NULL);
	int temperature = atoi(argv[4]); // Temperature
	
	int self_cutoff = 2;
	int bufferLength = 255;
	char buffer[bufferLength];
	char buffer2[bufferLength];
	
	double coordinates[nr_beads][3][nr_polymers];
	int contacts[((int) nr_atoms*(nr_atoms-1))/2][4];
	memset(coordinates,0.0,nr_atoms*3*sizeof(double));
	memset(contacts,0,2*((int) nr_atoms*(nr_atoms-1))*sizeof(int));
	
	
	snprintf(buffer2,sizeof(char)*bufferLength,"motifs_single_%dK_%.2f.txt",temperature,cutoff);
	//snprintf(buffer2,sizeof(char)*bufferLength,"motifs_%dK_%.2f.txt",temperature,cutoff);
	FILE *fp_motifs = fopen(buffer2,"w");
	if(fp_motifs == NULL ) {
    	perror ("Error opening file");
   		return(-1);
   	}
   	
   	snprintf(buffer2,sizeof(char)*bufferLength,"loops_single_%dK_%.2f.txt",temperature,cutoff);
   	//snprintf(buffer2,sizeof(char)*bufferLength,"loops_%dK_%.2f.txt",temperature,cutoff);
	FILE *fp_loops = fopen(buffer2,"w");
	if(fp_loops == NULL ) {
    	perror ("Error opening file");
   		return(-1);
   	}
   	
   	snprintf(buffer2,sizeof(char)*bufferLength,"contacts_single_%dK_%.2f.txt",temperature,cutoff);
   	//snprintf(buffer2,sizeof(char)*bufferLength,"contacts_%dK_%.2f.txt",temperature,cutoff);
	FILE *fp_contacts = fopen(buffer2,"w");
	if(fp_contacts == NULL ) {
    	perror ("Error opening file");
   		return(-1);
   	}
   	
   	
	for(int k=1; k<=100; k++){
		snprintf(buffer2,sizeof(char)*bufferLength,"./single_polymer/%dK%03d",temperature,k);
		//snprintf(buffer2,sizeof(char)*bufferLength,"./multiple_polymers/%dK%03d",temperature,k);
		FILE* fp_lammps = fopen(buffer2, "r");
			if(fp_lammps == NULL ){
		    		perror ("Error opening file");
		   		return(-1);
		   	}
		fgets(buffer, bufferLength, fp_lammps);
		
		while(fgets(buffer, bufferLength, fp_lammps)){
			char *token = strtok(buffer," ");
			char *ptr;
			int i = 0;
			int position = 0;
			int poly = 0;
			double cox = 0.0;
			double coy = 0.0;
			double coz = 0.0;
			
			while( token != NULL ){
				double strtoken = strtod(token,&ptr);
				//printf("%lf\t%d\n",strtoken,i);
				switch(i){
					case 0:
						position = (int) (((int) strtoken-1) % nr_beads);
						break;
					case 1: 
						poly = ((int) strtoken)-1;
						break;
					case 2:
						break;
					case 3:
						cox = strtoken;
						break;
					case 4:
						coy = strtoken;
						break;
					case 5: 
						coz = strtoken;
						break;
					default:
						perror ("Error while reading the line...");
						
				}			
	      			token = strtok(NULL, " ");
	      			i++;
	      			
	   		}
	   		coordinates[position][0][poly] = cox;
			coordinates[position][1][poly] = coy;
			coordinates[position][2][poly] = coz;
		}
		fclose(fp_lammps);
		
		double bl = bond_length(nr_beads, nr_polymers, coordinates);
	   	
		int poly1 = nr_polymers+1;
		int poly2 = nr_polymers+1;
		int pos1 = 0;
		int pos2 = 0;
		int count = 0;
		double xdist, ydist, zdist, dist;
		int intra=0, inter=0;
		
		for(int i=0; i<nr_atoms; i++){
			poly1 = floor(((double) i)/((double) nr_beads));
			pos1 = i-poly1*nr_beads;
			for(int j=i; j<nr_atoms; j++){
				poly2 = floor(((double) j)/((double) nr_beads));
				pos2 = j-poly2*nr_beads;
				xdist = (coordinates[pos1][0][poly1]-coordinates[pos2][0][poly2])*(coordinates[pos1][0][poly1]-coordinates[pos2][0][poly2]);
				ydist = (coordinates[pos1][1][poly1]-coordinates[pos2][1][poly2])*(coordinates[pos1][1][poly1]-coordinates[pos2][1][poly2]);
				zdist = (coordinates[pos1][2][poly1]-coordinates[pos2][2][poly2])*(coordinates[pos1][2][poly1]-coordinates[pos2][2][poly2]);
				dist = sqrt(xdist+ydist+zdist);
				
				if(poly1==poly2){ // Both on same polymer, take care of self-interaction cutoff!
					if((j-i>=self_cutoff)&&(dist<cutoff)){ // Intrachain contact is formed
						contacts[count][0] = poly1;
						contacts[count][1] = pos1;
						contacts[count][2] = poly2;
						contacts[count][3] = pos2;
						intra++;
						count++;
					}
				} else { // Both on different polymers, no problem
					if(dist<cutoff){ // Interchain contact is formed
						contacts[count][0] = poly1;
						contacts[count][1] = pos1;
						contacts[count][2] = poly2;
						contacts[count][3] = pos2;
						inter++;
						count++;
					}
				
				}
			}
		}
		int nr_contacts = count;
		
		
		int ct[11], ls[7] ={0, 0, 0, 0, 0, 0, 0}; // Loops only for (S P X I2 L2) (T2 I3)
		multichain_topology(bl, nr_contacts, nr_beads, nr_polymers, nr_atoms, contacts,coordinates, ct,ls);
		
		fprintf(fp_motifs,"%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",cutoff,ct[0],ct[1],ct[2],ct[3],ct[4],ct[5],ct[6],ct[7],ct[8],ct[9],ct[10],intra,inter);
		fprintf(fp_loops,"%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",cutoff,ls[0],ls[1],ls[2],ls[3],ls[4],ls[5],ls[6]);
		if(k==50){ // Use this specific simulation to export the contacts in order to build a contact map
			contact_map(fp_contacts, nr_contacts, nr_beads, nr_polymers, nr_atoms, contacts,coordinates, ct,ls);
		}	
	}
	fclose(fp_motifs);
	fclose(fp_loops);
	fclose(fp_contacts);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Execute time: %lfs \n", cpu_time_used); 
}


// Topology: S P X L2 T2 I2 T3 I3 I4 Unknown/Error
// Code:     1 2 3 4  5  6  7  8  9  0
void multichain_topology(double bl, int nr_contacts, int nr_beads, int nr_polymers, int nr_atoms, int contacts[((int) nr_atoms*(nr_atoms-1))/2][4], double coordinates[nr_beads][3][nr_polymers], int *ct, int *ls){
	int ct_motif = 0;
	int S=0, P=0, X=0, L2=0, T2=0, I2=0, T3=0, I3=0, I4=0, Err = 0;
	int poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4;
	
	for(int i=0; i<nr_contacts; i++){
		poly1 = contacts[i][0];
		pos1 = contacts[i][1];
		poly2 = contacts[i][2];
		pos2 = contacts[i][3];
		ct_motif = 0;
		
		for(int j=i+1; j<nr_contacts; j++){
			poly3 = contacts[j][0];
			pos3 = contacts[j][1];
			poly4 = contacts[j][2];
			pos4 = contacts[j][3];
			ct_motif= motif(poly1,poly2,poly3,poly4,pos1,pos2,pos3,pos4);
			
			switch(ct_motif){
				case 0:
					Err++;
					break;
				case 1:
					S++;
					ls[0]+=loop_size(nr_beads, nr_polymers, poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4, ct_motif);
					break;
				case 2:
					P++;
					ls[1]+=loop_size(nr_beads, nr_polymers, poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4, ct_motif);
					break;
				case 3:
					X++;
					ls[2]+=loop_size(nr_beads, nr_polymers, poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4, ct_motif);
					break;
				case 4:
					L2++;
					ls[4]+=loop_size(nr_beads, nr_polymers, poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4, ct_motif);
					break;
				case 5:
					T2++;
					ls[5]+=loop_size(nr_beads, nr_polymers, poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4, ct_motif);
					break;
				case 6:
					I2++;
					ls[3]+=loop_size(nr_beads, nr_polymers, poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4, ct_motif);
					break;
				case 7:
					T3++;
					break;
				case 8:
					I3++;
					ls[6]+=loop_size(nr_beads, nr_polymers, poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4, ct_motif);
					break;
				case 9:
					I4++;
					break;
				default:
					Err++;
					break;
			}
		}
	}
	ct[0] = Err; ct[1] = S; ct[2] = P; ct[3] = X; ct[4] = L2; ct[5] = T2; ct[6] = I2; ct[7] = T3; ct[8] = I3; ct[9] = I4; ct[10] = ((int) nr_contacts*(nr_contacts-1)/2);
}

void contact_map(FILE *fp, int nr_contacts, int nr_beads, int nr_polymers, int nr_atoms, int contacts[((int) nr_atoms*(nr_atoms-1))/2][4], double coordinates[nr_beads][3][nr_polymers], int *ct, int *ls){
	int poly1, poly2, poly3, poly4, pos1, pos2, pos3, pos4, ct_motif;
	
	for(int i=0; i<nr_contacts; i++){
		poly1 = contacts[i][0];
		pos1 = contacts[i][1];
		poly2 = contacts[i][2];
		pos2 = contacts[i][3];
		ct_motif = 0;
		for(int j=i+1; j<nr_contacts; j++){
			poly3 = contacts[j][0];
			pos3 = contacts[j][1];
			poly4 = contacts[j][2];
			pos4 = contacts[j][3];
			ct_motif= motif(poly1,poly2,poly3,poly4,pos1,pos2,pos3,pos4);
			fprintf(fp,"%d\t%d\t%d\n",i,j,ct_motif);
		}
	}
}

int motif(int poly1, int poly2, int poly3, int poly4, int pos1, int pos2, int pos3, int pos4){
int motif = 0;
	// First check first contact (1--2)
	if(poly1==poly2){ // Possible combinations: [S P X T2 I2 I3]
		if(poly3==poly4){ // Possible combinations: [S P X I2]
			if(poly3==poly1){ // Possible combinations: [S P X] // If all on the same polymer, pos1 is neccessarily smallest
				if(pos1<=pos3){ // Contact one is first
					if(pos2<=pos3){ // Resulting combination: [S]
						motif = 1;
					} else if((pos2>=pos3)&&(pos2<=pos4)){ // Resulting combination: [X]
						motif = 3;
					} else if((pos2>=pos3)&&(pos2>=pos4)){ // Resulting combination: [P]
						motif = 2;
					} else {
						motif = 0;
					}	
				} else if(pos3<=pos1){ // Contact two is first
					if(pos4<=pos1){ // Resulting combination: [S]
						motif = 1;
					} else if((pos4>=pos1)&&(pos4<=pos2)){ // Resulting combination: [X]
						motif = 3;
					} else if((pos4>=pos1)&&(pos4>=pos2)){ // Resulting combination: [P]
						motif = 2;
					} else {
						motif = 0;
					}
				} else {
					motif = 0;
				}
			} else { // Resulting combination: [I2]
				motif = 6;
			}
		} else { // Possible combinations: [T2 I3]
			if((poly3==poly1)||(poly4==poly1)){ // Resulting combinations: [T2]
				motif = 5;
			} else { // Resulting combinations: [I3]
				motif = 8;
			}
		}
	} else if(poly3==poly4){ // Possible combinations: [T2 I3]
		if((poly1==poly3)||(poly2==poly3)){// Resulting combinations: [T2]
			motif = 5;
		} else { // Resulting combinations: [I3]
			motif = 8;
		}
	
	} else { // Possible combinations: [L2 T3 I4]
		if((poly1!=poly3)&&(poly1!=poly4)&&(poly2!=poly3)&&(poly2!=poly4)){ // Resulting combinations: [I4]
			motif = 9;
		} else if(((poly1==poly3)&&(poly2==poly4))||((poly1==poly4)&&(poly2==poly3))){ // Resulting combinations: [L2]
			motif = 4;
		} else { // Resulting combinations: [T3]
			motif = 7;
		}
	}
	return motif;
}

int loop_size(int nr_beads, int nr_polymers, int poly1, int poly2, int poly3, int poly4, int pos1, int pos2, int pos3, int pos4, int ct_motif){
	int L;
	switch(ct_motif){
		case 0:
			break;
		case 1: // S motif; size is determined by individual loops on same polymer
			L = (abs(pos2-pos1) + abs(pos4-pos3));
			break;
		case 2: // P motif; size is determined by two pinched loops but total loop size is just from smallest position to largest
			if(pos1<pos3){ // Order is 1-34-2 -- ABBA
				L = abs(pos2-pos1);
			} else { // Order is 3-12-4 -- BAAB
				L = abs(pos4-pos3);
			}
			break;
		case 3: // X motif; size is determined by two pinched loops but total loop size is just from smallest position to largest
			if(pos1<pos3){ // Order is 1324 -- ABAB
				L = abs(pos4-pos1);
			} else { // Order is 3142 -- BABA
				L = abs(pos2-pos3);
			
			}
			break;
		case 4: // L2 motif; size is determined by 'bubble' between two chains
			if((poly1==poly3)&&(poly2==poly4)){ // 1 and 3 share a chain, and 2 and 4 another
				L = abs(pos1-pos3) + abs(pos2-pos4);
			}
			break;
		case 5: // T2 motif; only single loop determines size
			if(poly1==poly2){ // Loop is between contacts 1-2
				L = abs(pos2-pos1);
			} else if(poly3==poly4){
				L = abs(pos4-pos3);
			} else {
				L=0;
			}
			break;
		case 6: // I2 motif; two independent loops
			L = abs(pos2-pos1) + abs(pos4-pos3);
			break;
		case 7: // T3 motif; no loops are formed 
			break;
		case 8: // I3 motif; again only single loop determines size
			if(poly1==poly2){ // Loop is on poly1
				L = abs(pos2-pos1);
			} else if(poly3==poly4){ // Loop is on poly3
				L = abs(pos4-pos3);
			} else {
				L=0;
			}
			break;
		case 9: //I4 motif; no loops
			break;
		default:
			break;
	}
	return L;
}

double bond_length(int nr_beads, int nr_polymers, double coordinates[nr_beads][3][nr_polymers]){
	double xdist = (coordinates[0][0][0]-coordinates[1][0][0])*(coordinates[0][0][0]-coordinates[1][0][0]);
	double ydist = (coordinates[0][1][0]-coordinates[1][1][0])*(coordinates[0][1][0]-coordinates[1][1][0]);
	double zdist = (coordinates[0][2][0]-coordinates[1][2][0])*(coordinates[0][2][0]-coordinates[1][2][0]);
	return sqrt(xdist+ydist+zdist);
}
