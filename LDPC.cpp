#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <windows.h>
#include <math.h>
#include <algorithm>

#define debug 0

int frame = 0;
int frame2 = 0;
const int length = 10;
const int checker = 5;
const int error_target = 200;
int times;
double stddev;
double Eb_N0_dB;
double code_rate = 1./1.;

int debug_message[length] = {-1,1,1,-1,1,-1,1,1,-1,1};

using namespace::std;

double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

int sign(double a)
{
    return a>=0.?1:-1;
}

const int P[checker][length] = 
{
	{1,1,0,0,1,1,0,0,0,0},
	{0,0,1,1,0,1,1,0,0,0},
	{1,0,0,0,0,0,1,1,1,0},
	{0,1,1,0,0,0,0,1,0,1},
	{0,0,0,1,1,0,0,0,1,1}
};

struct bit
{
	int value;
	double prob;
};

void change(bit &a)
{
	if (a.value == 1)	a.value = -1;
	else				a.value = 1;
	a.prob = 1-a.prob;
}

void prob_calc(bit* result, double* array, double stddev)
{
	for (int i=0;i<length;i++)
	{
		if (array[i]>=0)
		{
			result[i].value = 1;
			result[i].prob  = exp(-0.5*(array[i]-1.)/stddev*(array[i]-1.)/stddev)/(exp(-0.5*(array[i]-1.)/stddev*(array[i]-1.)/stddev)+exp(-0.5*(array[i]+1.)/stddev*(array[i]+1.)/stddev));
		}
		else
		{
			result[i].value = -1;
			result[i].prob  = exp(-0.5*(array[i]+1.)/stddev*(array[i]+1.)/stddev)/(exp(-0.5*(array[i]-1.)/stddev*(array[i]-1.)/stddev)+exp(-0.5*(array[i]+1.)/stddev*(array[i]+1.)/stddev));
		}
		if (debug) printf ("array %+lf, value %+d, prob %lf\n",array[i],result[i].value,result[i].prob);
	}
}

void prob_calc(bit* result, int* array, double stddev)
{
	for (int i=0;i<length;i++)
	{
		if (array[i]>=0)
		{
			result[i].value = 1;
			result[i].prob  = 1;
		}
		else
		{
			result[i].value = -1;
			result[i].prob  = 1;
		}
		if (debug) printf ("array %+d, value %+d, prob %lf\n",array[i],result[i].value,result[i].prob);
	}
}

bool parity_check(bit* array, int* result)
{
	if (debug)	printf ("Parity check: ");
	bool r = true;
	int notpass;
	for (int i=0;i<checker;i++)
	{
		notpass = 0;
		for (int j=0;j<length;j++)
		{
			if (P[i][j])
			{
				notpass+=array[j].value==1 ? 1 : 0;
			}
		}
		if (notpass%2)
		{
			result[i] = 1;
			r = false;
		}	
		else			result[i] = 0;
		if (debug)	printf ("%d ",result[i]);
	}
	if (debug)	printf ("\n");
	return r;
}

void LDPC(bit* _array, int* _visited,double my_rate,double* rate,int* success_codeword)
{
	if(my_rate < *rate)	return;
	bit* array;
	array = new bit [length];
	int* visited;
	visited = new int [length];
	for (int i=0;i<length;i++)
	{
		array[i] = _array[i];
		visited[i] = _visited[i];
	}
	int* notpass;
	notpass = new int [checker];
	
	parity_check(array, notpass);
	
	bool success_flag = true;
	for (int i=0;i<checker;i++)
	{
		if (notpass[i])
		{
			success_flag = false;
			break;
		}	
	}
	if (success_flag)
	{
		if (debug) printf ("Success codeword\n");
		for (int i=0;i<length;i++)
		{
			if (debug) printf ("%d ",array[i].value);
		}
		if (debug) printf ("\n");
		if (debug) printf ("rate = %.20lf\n",my_rate);
		if (my_rate > *rate)
		{
			*rate = my_rate;
			for (int i=0;i<length;i++)
			{
				success_codeword[i]=array[i].value;
			}
		}
		
		//char pause;
		//scanf (" %c",&pause);
	}
	
	for (int i=0;i<checker;i++)
	{
		if (notpass[i])
		{
			for (int j=0;j<length;j++)
			{
				if (P[i][j] && visited[j]==0)
				{
					visited[j] = 1;
					change(array[j]);
					if (debug)	printf ("visit %d \n",j);
					LDPC(array,visited,my_rate*array[j].prob,rate,success_codeword);
					change(array[j]);
					visited[j] = 0;
				}
			}
		}
	}
	
	delete [] notpass;
	delete [] visited;
	delete [] array;
}

double avg ( vector <double> & v )
{
    double return_value = 0.0;
    int n = v.size();
       
    for (int i=0;i<n;i++)
    {
    	return_value += v[i];
    }
       
    return (return_value / n);
}
//****************End of average funtion****************


//Function for variance
double variance (vector <double> & v , double mean)
{ 
    double sum = 0.0;
    double var =0.0;
    int n = v.size();
       
    for (int j=0;j<n;j++)
    {
        sum += pow((v[j]-mean),2);
    }
    var = sum/(n-1);
    return var;
}
//****************End of variance funtion****************

int main()
{
    srand(time(NULL));
    int* message_array;
    double* noise_array;
    double* receive_array;
    
    vector <double> result_rate;
    vector <double> LDPC_time;

    printf ("length = %d\n",length);

	char filename[1024];
	char filename_hard_decision[1024];
	char filename_rate[1024];
	char filename_time[1024];
	sprintf (filename,"LDPC_%d.txt",length);
	FILE *f = fopen(filename,"w");
	fclose(f);
	sprintf (filename_hard_decision,"Hard_Decision_%d.txt",length);
	f = fopen(filename_hard_decision,"w");
	fclose(f);
	sprintf (filename_rate,"Rate_%d.txt",length);
	f = fopen(filename_rate,"w");
	fclose(f);
	sprintf (filename_time,"time_%d.txt",length);
	f = fopen(filename_time,"w");
	fclose(f);

	for (Eb_N0_dB=-3.;Eb_N0_dB<6.41;Eb_N0_dB+=0.1)
	{
		stddev = sqrt(pow(10,-Eb_N0_dB/10)/2/code_rate);
		times = 1;																				// change here
		cout <<"SNR = "<<Eb_N0_dB<<", stddev = "<<stddev<<endl;
		while (times)
		{
		    int error_count = 0;
		    int error_count2 = 0;
		    //int error_count_bit = 0;
		    frame = 0;
		    frame2 = 0;
		    while(1)
		    {
			    message_array = new int[length];
			    noise_array = new double[length];
			    receive_array = new double[length];
			    
			    bit* array;
			    array = new bit [length];
			    int* visited;
			    visited = new int [length];
			    
			    double my_rate = 1;
				double rate = 0;
				int* success_codeword;
				success_codeword = new int [length];

			    for (int i=0;i<length;i++)      //construct initial message
			    {
			        message_array[i] = rand()/0.5 > RAND_MAX ? 1 : -1;
			        if (debug)	message_array[i] = debug_message[i];
			    }
			    
			    for (int i=0;i<length;i++)      
			    {
			        if (debug) printf ("%d ",message_array[i]);
			    }
			    if (debug) printf ("\n");
			    
			    prob_calc(array,message_array,stddev);
			    int* notpass;
				notpass = new int [checker];
				bool success_generate = parity_check(array, notpass);
				if (!success_generate)
				{
					goto enddelete;
				}	

				frame++;
		    	frame2++;

				if (debug)	printf ("codeword generate, start to add noise\n");
			
			    for (int i=0;i<length;i++)
			    {
			        noise_array[i] = rand_normal(0,stddev);
			        receive_array[i] = message_array[i] + noise_array[i];
			        if (debug) cout << receive_array[i] << " ";
			    }
			    if (debug) printf ("\n");

				prob_calc(array,receive_array,stddev);
				for (int i=0;i<length;i++)
					visited[i] = 0;
				
				LARGE_INTEGER t1, t2, ts;
			    QueryPerformanceFrequency(&ts);
			    QueryPerformanceCounter(&t1);
			    
				LDPC(array,visited,my_rate,&rate,success_codeword);
				
				
			    QueryPerformanceCounter(&t2);
			    if (rate<0.99)
			    {
			    	result_rate.push_back(rate);
					LDPC_time.push_back(1000000.*(t2.QuadPart-t1.QuadPart)/(double)(ts.QuadPart));
				}
				
				
				//printf ("%lf %lf\n",rate,1000000.*(t2.QuadPart-t1.QuadPart)/(double)(ts.QuadPart));
				//char pause;
				//scanf ("%c",&pause);
				
				//Summary
			    for (int i=0;i<length;i++)
			    {
			        if (success_codeword[i]!=message_array[i])
					{
						error_count++;
						break;					//Frame error rate
					}
				}
				for (int i=0;i<length;i++)
			    {
			        if (sign(receive_array[i])!=message_array[i])
					{
						error_count2++;
						break;					//Frame error rate
					}
				}
				
			enddelete:
				
			    delete [] message_array;
			    delete [] receive_array;
			    delete [] notpass;
			    delete [] array;
			    delete [] visited;
			    delete [] noise_array;
			    delete [] success_codeword;
			    
			    if (error_count>=error_target)	break;
			}
			printf ("frame count = %d, ",frame);

			double avg_rate = avg(result_rate);
			double stddev_rate = sqrt(variance(result_rate,avg_rate));
			
			f = fopen(filename_rate,"a");
			fprintf(f,"%.1lf %lf %lf",(double) Eb_N0_dB,avg_rate,stddev_rate);
			fprintf(f,"\n");
			fclose(f);
			result_rate.clear();
			
			sort(LDPC_time.begin(),LDPC_time.end());
			double avg_time = avg(LDPC_time);
			double stddev_time = sqrt(variance(LDPC_time,avg_time));
			
			f = fopen(filename_time,"a");
			fprintf(f,"%.1lf %lf %lf %lf %lf",(double) Eb_N0_dB,avg_time,stddev_time,LDPC_time.front(),LDPC_time.back());
			fprintf(f,"\n");
			fclose(f);
			LDPC_time.clear();

			f = fopen(filename,"a");
			fprintf (f,"%.1lf %d %d\n",(double) Eb_N0_dB,error_count,frame);					//frame error rate
			fclose(f);
			
			printf ("frame count = %d\n",error_count2);

			f = fopen(filename_hard_decision,"a");
			fprintf (f,"%.1lf %d %d\n",(double) Eb_N0_dB,error_count2,frame2);					//frame error rate
			fclose(f);

			char temp;
			if (debug) scanf ("%c",&temp);

			times--;
		}
	}


    return 0;
}

