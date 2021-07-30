#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stddef.h>
#include <mpi.h>
#define NBR_CITIES 52
#define POP_SIZE 200
#define NBR_GEN 200


char * DataSetFileName;

typedef struct city{
        int id;
        int x;
        int y;
}city;

typedef struct tour{
        city cities[NBR_CITIES];
        float fitness;
}tour;

typedef struct population{
        tour tours[POP_SIZE];
        tour fittest;
}population;

// List of Global Varibales
static tour OriginalTour;
static int elitismSizeG, touranmentPoolSizeG;
static float mutationRateG;
static int Population_size;
static float SELECTION_PROBABILITY = 0.1;
static int generationCounter = 0;

//// List Of All Methods
// DataSet Methods
void readDataSet(char* FileName);
tour getTour();
// Tour Methods
float calculate_Fitness(tour t1);
int containCity(tour t1,city c);
int compare(tour t1,tour t2);
tour getRandomTour(tour t);
void toString(tour t);
int getIndexOfCity(tour t , city c);
// Population Methods
void init_population(population *pop);
tour getFittestTour(population pop);
tour* getNFittestTours(population pop , int n);
void sortPopulationASC(population *pop);
void sortPopulationDESC(population *pop);
void toStringPopulation(population pop);
tour* sortToursASC(tour t[] ,int size);
void calculateFitnessForAll(population *p);

//Function to define Our Types
void defineTypes();

// Geenetic Algorithm Methods
void GeneticAlgorithm(int numberOfGenerations,int elitismSize , float mutationRate , int touranmentPoolSize,int argc,char** argv);
population reproduction(population p);
//Mutation Methods
void swapMutation(tour *tour);
//CrossOverMethods
tour OnePointCrossOver(tour parent1 , tour parent2);
tour TwoPointCrossOver(tour parent1 , tour parent2);
tour ox1CrossOver(tour parent1 , tour parent2);
//Selection Methods
tour* elitismSelection(population p);
tour touranmentSeelection(population p);
tour rouletteWheelSeelction(population p);
tour RankSeelction(population p);
//function to Write The Result in a file
void WriteResultIntoFile(population pop);

// Function To Read The DataSet
void readDataSet(char* FileName)
{
        FILE *Pfile;
        int i=0;

        Pfile = fopen(FileName,"r");

        if(Pfile == NULL)
        {
                printf("Failed To Open The File . \n");
                return;
        }

        while(!feof(Pfile))
        {
                fscanf(Pfile,"%d\t%d\t%d\n",&OriginalTour.cities[i].id,&OriginalTour.cities[i].x,&OriginalTour.cities[i].y);
                i++;
        }
        OriginalTour.fitness = calculate_Fitness(OriginalTour);
        fclose(Pfile);
}

tour getTour()
{
                return OriginalTour;
}

//Fitness Function
float calculate_Fitness(tour t1)
{
        float fitness = 0;
        int i;
        for(i=0;i<NBR_CITIES;i++)
        {
                fitness += sqrt(pow(t1.cities[i+1].x - t1.cities[i].x, 2) + pow(t1.cities[i+1].y - t1.cities[i].y, 2)) ;
        }
        fitness += sqrt(pow(t1.cities[NBR_CITIES - 1].x - t1.cities[0].x, 2) + pow(t1.cities[NBR_CITIES - 1].y - t1.cities[0].y, 2)) ;
        return fitness;
}
//Function to Verify if a Tour Contain a City
int containCity(tour t1,city c)
{
        int i,contain = 0;
        for(i=0;i<NBR_CITIES;i++)
        {
                if(t1.cities[i].id == c.id)
                        contain = 1;
        }
        return contain;
}
//Function to Compare two tours
int compare(tour t1,tour t2)
{
        return t1.fitness > t2.fitness ? 1 : t1.fitness < t2.fitness ? -1 : 0 ;
}
//Function to get a random Tour for a tour passed in The Parameters
tour getRandomTour(tour t)
{
        int index;
                for (int i = 1; i < NBR_CITIES - 1; i++)
        {
          index = (int)rand()%(-50) + 2;
                  city c = t.cities[index];
          t.cities[index] = t.cities[i];
          t.cities[i] = c;
        }
        t.fitness = calculate_Fitness(t);
        return t;
}
//Function to Display a Tour
void toString(tour t)
{
        for(int i = 0 ; i < NBR_CITIES ; i++)
        {
                printf("%d",t.cities[i].id);
                if(i != NBR_CITIES -1)
                {
                        printf("->");
                }
        }
        printf("\n");
}

//Function to get The Index of a specific City
int getIndexOfCity(tour t , city c)
{
        for(int i=0;i<NBR_CITIES ; i++)
        {
                if(t.cities[i].id == c.id)
                {
                        return i;
                }
        }
}
//Function to initialize The Population
void init_population(population *pop)
{
        readDataSet(DataSetFileName);
        tour tour = getTour(),*pointerT;
        pop->tours[0] = tour;
        for(pointerT = pop->tours + 1 ; pointerT < pop->tours + POP_SIZE ; pointerT++)
        {
                *(pointerT) = getRandomTour(*(pointerT - 1));
        }
        calculateFitnessForAll(pop);
        pop->fittest = getFittestTour(*pop);
}
//Function To Get The Fittest Tour In a Population
tour getFittestTour(population pop)
{
        sortPopulationASC(&pop);
        return pop.tours[0];
}
//Function to Get Top N Fittest Chromosome in a Population
tour* getNFittestTours(population pop,int n)
{
        tour *array = (tour*) malloc(n * sizeof(tour));
        sortPopulationASC(&pop);
        for(int i=0;i<n;i++)
        {
                array[i] = pop.tours[i];
        }
        return array;
}
//Function to Sort The Population ASC based on The Fitness of their Tours
void sortPopulationASC(population *pop)
{
        calculateFitnessForAll(pop);
        for(int i=0;i<POP_SIZE-1 ; i++)
        {
                for(int j=i+1;j<POP_SIZE;j++)
                {
                        if(compare(pop->tours[i],pop->tours[j]) > 0)
                        {
                                tour tmp = pop->tours[i];
                                pop->tours[i] = pop->tours[j];
                                pop->tours[j] = tmp;
                        }
                }
        }
}
//Function to Sort The Population DESC based on The Fitness of their Tours
void sortPopulationDESC(population *pop)
{
        for(int i=0;i<POP_SIZE-1 ; i++)
        {
                for(int j=i+1;j<POP_SIZE;j++)
                {
                        if(compare(pop->tours[i],pop->tours[j]) < 0)
                        {
                                tour tmp = pop->tours[i];
                                pop->tours[i] = pop->tours[j];
                                pop->tours[j] = tmp;
                        }
                }
        }
}
//Function To Display The Tours of a population
void toStringPopulation(population pop)
{
        for(int i=0;i<POP_SIZE ; i++)
        {
                printf("Tour %d ,  Fitness = %.2f \n",i+1,calculate_Fitness(pop.tours[i]));
        }
}

//Function to sort a set of tours ASC
tour* sortToursASC(tour t[] ,int size)
{
        for(int i=0;i<size-1 ; i++)
        {
                for(int j=i+1;j<size;j++)
                {
                        if(compare(t[i],t[j]) == 1)
                        {
                                tour tmp = t[i];
                                t[i] = t[j];
                                t[j] = tmp;
                        }
                }
        }
        return &t[0];
}
//Function to calculate The Fitness for all the elements of a Population
void calculateFitnessForAll(population *p)
{
        for(int i=0;i<POP_SIZE ; i++)
        {
                p->tours[i].fitness = calculate_Fitness(p->tours[i]);
        }
}

void myReduceFunction(void *in , void *out , int *len,MPI_Datatype *type)
{
        tour *invalue = (tour *)in;
        tour *outvalue = (tour *)out;

        outvalue[0] = invalue[0] ;


        for(int i=1;i<*len;i++)
        {
                if(invalue[i].fitness < outvalue[0].fitness)
                {
                        outvalue[0] = invalue[i];
                }
        }
        return;
}

void swap (tour *a, tour *b)
{
        tour temp = *a;
        *a = *b;
        *b = temp;
}

void shuffle (tour arr[], int n)
{
        srand (time(NULL));
        for (int i = n - 1; i > 0; i--)
        {
                int j = rand() % (i + 1);
                swap(&arr[i], &arr[j]);
        }
}

// Genetic Algorithm
void GeneticAlgorithm(int numberOfGenerations,int elitismSize , float mutationRate , int touranmentPoolSize,int argc,char** argv)
{
        elitismSizeG = elitismSize;
        mutationRateG = mutationRate ;
        touranmentPoolSizeG = touranmentPoolSize ;


        int procID , nbrProcs,code,next,prev,sendSize;
        const int tag = 10;
        MPI_Status status;
        tour tourstoSend[POP_SIZE - elitismSizeG];
        //Initialize The Parallel Region
        MPI_Init(&argc,&argv);
        MPI_Comm_size(MPI_COMM_WORLD,&nbrProcs);
        MPI_Comm_rank(MPI_COMM_WORLD,&procID);


        //Define Our Types
        //define City Type
        const int nitems1=3;
        int blocklengths1[3] = {1,1,1};
        MPI_Datatype types1[3] = {MPI_INT, MPI_INT,MPI_INT};
        MPI_Datatype city_type;
        MPI_Aint     offsets1[3];

        offsets1[0] = offsetof(city, id);
        offsets1[1] = offsetof(city, x);
        offsets1[2] = offsetof(city, y);

        MPI_Type_create_struct(nitems1, blocklengths1, offsets1, types1, &city_type);
        MPI_Type_commit(&city_type);
        //define Tour Type
        const int nitems2=2;
        int blocklengths2[2] = {NBR_CITIES,1};
        MPI_Datatype types2[2] = {city_type, MPI_FLOAT};
        MPI_Datatype tour_type;
        MPI_Aint     offsets2[2];
        offsets2[0] = offsetof(tour, cities);
        offsets2[1] = offsetof(tour, fitness);
        MPI_Type_create_struct(nitems2, blocklengths2, offsets2, types2, &tour_type);
        MPI_Type_commit(&tour_type);
        //End of types defintion

        //Define MPI OP
        MPI_Op myReduceOperation;
        MPI_Op_create(myReduceFunction,1,&myReduceOperation);
        //End Definition

        //Open The File to Write The Final Result
        FILE *pf;
        pf = fopen("parallelData2.txt","a+");
        if(pf == NULL)
        {
                printf("Failed To Open The File to Write The Result .\n");
                return ;
        }

        //Initialize The Population
        population initial_Population ;

        init_population(&initial_Population);
        tour TheBestTour=initial_Population.tours[0],tmp;
        TheBestTour.fitness = calculate_Fitness(TheBestTour);
        if(procID == 0)
        {
        fprintf(pf,"%d\t%f\n",generationCounter,TheBestTour.fitness);
        printf("The Initial Population , Fitness : %f \n",getFittestTour(initial_Population).fitness);
        toString(getFittestTour(initial_Population));
        }
        float epidemic_prob;


        population newPopulation = initial_Population;
        for(int i=0;i<numberOfGenerations ; i++)
        {
                newPopulation = reproduction(newPopulation);
                calculateFitnessForAll(&newPopulation);
                generationCounter ++;
                mutationRateG += 0.0009;
                //printf("Generation Number : %d , Fitness : %f \n",generationCounter,newPopulation.fittest.fitness);
                sortPopulationASC(&newPopulation);
                for(int i=elitismSizeG,j=0;i<POP_SIZE;i++,j++)
                {
                        tourstoSend[j] = newPopulation.tours[i];
                }
                sendSize = (POP_SIZE - elitismSizeG);
                shuffle(tourstoSend,sendSize);
                next = (procID + 1)%nbrProcs ;
                prev = (nbrProcs + procID -1)%nbrProcs ;

                if(procID == 0)
                {
                        MPI_Send(tourstoSend,sendSize,tour_type,next,tag,MPI_COMM_WORLD);

                        MPI_Recv(tourstoSend,sendSize,tour_type,prev,tag,MPI_COMM_WORLD,&status);
                }else{
                        MPI_Recv(tourstoSend,sendSize,tour_type,prev,tag,MPI_COMM_WORLD,&status);

                        MPI_Send(tourstoSend,sendSize,tour_type,next,tag,MPI_COMM_WORLD);
                }

                for(int i=elitismSizeG,j=0;i<POP_SIZE;i++,j++)
                {
                        newPopulation.tours[i] = tourstoSend[j];
                }
                calculateFitnessForAll(&newPopulation);
                MPI_Reduce(&newPopulation.fittest,&tmp,1,tour_type,myReduceOperation,0,MPI_COMM_WORLD);
                if(procID == 0)
                {
                        fprintf(pf,"%d\t%f\n",generationCounter,newPopulation.fittest.fitness);
                        tmp.fitness =   calculate_Fitness(tmp);
                        TheBestTour = tmp.fitness > TheBestTour.fitness ? TheBestTour : tmp;

                        newPopulation.tours[POP_SIZE -1] = TheBestTour;
                        calculateFitnessForAll(&newPopulation);
                        printf("Generation Number : %d , Fitness : %f \n",generationCounter,newPopulation.fittest.fitness);
                }
        }
        if(procID == 0)
        {
        printf("After %d Generation , The Final Result is : \n",generationCounter);
        printf("Fitness :  %f .\n",newPopulation.fittest.fitness);
        printf("The Final Tour is : \n");
        toString(newPopulation.fittest);
        WriteResultIntoFile(newPopulation);

        }
        fclose(pf);
        MPI_Finalize();


}

// Method to Write The Result in a File
void WriteResultIntoFile(population pop)
{
        FILE *pf;
        pf = fopen("parallelResult2.txt","w");
        if(pf == NULL)
        {
                printf("Failed to Open The File .\n");
        }

        for(int i=0;i<NBR_CITIES;i++)
        {
                fprintf(pf,"%d\t%d\t%d\n",pop.fittest.cities[i].id,pop.fittest.cities[i].x,pop.fittest.cities[i].y);
        }
        fclose(pf);
}

//Reproduction Method that Contain {Selection -> CrossOver -> Mutation}
population reproduction(population p)
{
        population newPopulation;
        int index;
        tour *elitism = elitismSelection(p);
        for(index = 0 ; index < elitismSizeG ; index ++)
        {
                newPopulation.tours[index] = *(elitism+index);
        }

        tour parent1,parent2,child,child1,child2;
        int selectionMethod,CrossOverMethod;
        for(int i=index;i<POP_SIZE;i++)
        {

                parent1 = touranmentSeelection(p);
                parent2 = touranmentSeelection(p);

                //One Point CrossOver
                child1 = OnePointCrossOver(parent1,parent2);
                child2 = OnePointCrossOver(parent1,parent2);

                child1.fitness = calculate_Fitness(child1);
                child2.fitness = calculate_Fitness(child2);

                child = child1.fitness > child2.fitness ? child2 : child1;
                //Two Point CrossOver
                child1 = TwoPointCrossOver(parent1,parent2);
                child2 = TwoPointCrossOver(parent1,parent2);

                child1.fitness = calculate_Fitness(child1);
                child2.fitness = calculate_Fitness(child2);

                child = child1.fitness > child2.fitness ? child2 : child1;
                //OX1 CrossOver
                child1 = ox1CrossOver(parent1,parent2);
                child2 = ox1CrossOver(parent1,parent2);

                child1.fitness = calculate_Fitness(child1);
                child2.fitness = calculate_Fitness(child2);

                child = child1.fitness > child2.fitness ? child2 : child1;


                child1 = child;
                swapMutation(&child1);
                child1.fitness = calculate_Fitness(child1);

                child = child.fitness > child1.fitness ? child1 : child;

                newPopulation.tours[i] = child;
        }
        calculateFitnessForAll(&newPopulation);
        newPopulation.fittest = getFittestTour(newPopulation);
        return newPopulation;
}

// Swap Method for Mutation
void swapMutation(tour *tour)
{
        double p;
        int mutationPoint ;

        for(int i=1;i<NBR_CITIES;i++)
        {
                        p = (double)rand() / (double) RAND_MAX ;
                        if(p < mutationRateG)
                        {
                                mutationPoint = (int)rand()%(-50) + 2;
                                city tmp = tour->cities[i];
                                tour->cities[i] = tour->cities[mutationPoint];
                                tour->cities[mutationPoint] = tmp;
                        }
        }
        tour->fitness = calculate_Fitness(*tour);
}

// Selection Methods
// Elitism Selection
tour* elitismSelection(population p)
{
        sortPopulationASC(&p);
        tour *pt,*t = (tour *) malloc(elitismSizeG * sizeof(tour));
        for(pt = t ; pt < t + elitismSizeG ; pt++)
        {
                *pt = p.tours[pt-t];
        }
        return t;
}

// Touranment Selection
tour touranmentSeelection(population p)
{
                        tour touranmentPool[touranmentPoolSizeG];
                        float Selection_prob;
                        int touranmentIndex = 0,index;
                        while(touranmentIndex < touranmentPoolSizeG)
                        {
                                index = (int) rand() % (-POP_SIZE) + 1;
                                Selection_prob = (float) rand() / (float) RAND_MAX;
                                if(Selection_prob < SELECTION_PROBABILITY)
                                {
                                        touranmentPool[touranmentIndex] = p.tours[index];
                                        touranmentIndex++;
                                }
                        }
                return *sortToursASC( touranmentPool, touranmentPoolSizeG);
}

//RouletteWheel Selection
tour rouletteWheelSeelction(population p)
{
                float sumOfFitnesses = 0;
                tour selectedTour;
                calculateFitnessForAll(&p);
                for(int i=0;i<POP_SIZE;i++)
                {
                        sumOfFitnesses += p.tours[i].fitness;
                }

                float count = (float) (((float)rand() / (float) RAND_MAX) * sumOfFitnesses);


                for(int i=0;i<POP_SIZE;i++)
                {
                        count += p.tours[i].fitness;
                        if(count > sumOfFitnesses)
                        {
                                selectedTour = p.tours[i];
                                break;
                        }
                }
                return selectedTour;
}

//Rank Selection
tour RankSeelction(population p)
{
                tour selectedTour ;
                sortPopulationASC(&p);

                int count = (int) (((float)rand() / (float) RAND_MAX) * POP_SIZE);

                for(int i=0;i<POP_SIZE;i++)
                {
                        count += i;
                        if(count > POP_SIZE)
                        {
                                selectedTour = p.tours[i];
                                break;
                        }
                }
                return selectedTour;
}

// CrossOver Methods
// One Point CrossOver
tour OnePointCrossOver(tour parent1 , tour parent2)
{
        int index,crossOverPoint = (int)rand()%(-50) + 2;

        for(int i=crossOverPoint;i<NBR_CITIES ; i++)
        {
                city c = parent2.cities[i];
                index = getIndexOfCity(parent1,c);
                parent1.cities[index] = parent1.cities[i];
                parent1.cities[i] = c ;
        }

        parent1.fitness = calculate_Fitness(parent2);

        return parent1;
}

// Two Point CrossOver
tour TwoPointCrossOver(tour parent1 , tour parent2)
{
        int index,crossOverPointStart,crossOverPointEnd;
        //Initialize The two CrossOver Points
        crossOverPointEnd = (int)rand()%(-50) + 2 ;
        crossOverPointStart = (int)rand()%(-50) + 2;
        //Make sure that the start index less than The End Index
        while(crossOverPointStart < crossOverPointEnd )
        {
                crossOverPointEnd = (int)rand()%(-50) + 2 ;
                crossOverPointStart = (int)rand()%(-50) + 2;
        }
        //Copy The elements of parent2 in parent1
        for(int i=crossOverPointStart;i<crossOverPointEnd ; i++)
        {
                city c = parent2.cities[i];
                index = getIndexOfCity(parent1,c);
                parent1.cities[index] = parent1.cities[i];
                parent1.cities[i] = c ;
        }

        parent1.fitness = calculate_Fitness(parent2);
        return parent1;
}

//OX1 CrossOver which is cycle crossover or Davis Order CrossOver
tour ox1CrossOver(tour parent1 , tour parent2)
{
        int Start,End;
        tour child;
        //Insert City 1 in Index 0
        child.cities[0] = parent1.cities[0];
        //Initialze The Tour With NULL
        city nullCity;
        nullCity.id = -1;
        for(int i=1;i<NBR_CITIES;i++)
        {
                child.cities[i] = nullCity;
        }

        //Calculate The Start and End indexes
        End = (int)rand()%(-50) + 2 ;
        Start = (int)rand()%(-50) + 2;
        //Make sure That Start Index less than End Index
        while(End <= Start)
        {
                End = (int)rand()%(-50) + 2 ;
                Start = (int)rand()%(-50) + 2;
        }

        //copy The element of Parent1 between Start and end Indexes in The Child
        for(int i=Start;i<End;i++)
        {
                        child.cities[i] = parent1.cities[i];
        }
        //Fill The Rest of Child with Cities in Parent2
        for(int i=1;i<NBR_CITIES;i++)
        {
                if(!containCity(child,parent2.cities[i]))
                {
                        for(int j=1;j<NBR_CITIES;j++)
                        {
                                if(child.cities[j].id == -1)
                                {
                                        child.cities[j] = parent2.cities[i];
                                        break;
                                }
                        }
                }
        }
        return child;
}



// ------------------------------------------ Main Method --------------------------------------------
int main(int argc , char** argv)
{
        if(argc < 3)
        {
                printf("Not Enough Parameters .\n");
        }
        DataSetFileName = argv[1];
        generationCounter = NBR_GEN * atoi(argv[2]);
        srand(time(NULL));
        GeneticAlgorithm(NBR_GEN,20,0.1,10,argc,argv);



        return 0;
}
