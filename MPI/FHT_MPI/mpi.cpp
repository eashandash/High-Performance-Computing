   #include <mpi.h>
   #include <iostream>
   #include <fstream>
   #include <sstream>
   #include <math.h> 
   #include <vector>
   #include <time.h>
   #include <cstdlib>
   #include <bitset> // we require it for special binary representation
  
  using namespace std;
  // input data 
  double * data_in_r;
  double * data_in_i;
  
  //output data
  double * data_out_r;
  double * data_out_i;
  
  //working buffer
  double * buffer_r;
  double * buffer_i;
  
  //reverse binary input data
  double * rb_buffer_r;
  double * rb_buffer_i;
 
  int p_length=240;
  int p2_length=0;
  int bitlength;
  int what_power=0;
  unsigned int mask=0x00000001;
  
 // MPI variables
 int world_rank;  // the rank of this process
 int world_size; // size of world, thus number of processors in this world
 int FROM;
 int TO;
 int rFROM;
 int rTO;

//First, we have to include again the mpi.h include file to be able to use the MPI environment. For the communication, we need some variables. As before, the world_rank and world_size variables will contain the number of the actual working node and the number of all nodes. For data partition, we define some variables, with which we can follow up the data range for different nodes (FROM, TO, rFROM, rTO).

//Since the parallel solution will use the same functions as the serial one, we define them the same way.
 
 unsigned int revbit(unsigned int val)
 {
 	unsigned int rval=0;
 	unsigned int fval=0;
 	for(int i=0;i<what_power;i++)
 	{
 	   rval = (val >> i) & mask;
 	   fval = fval | (rval << (what_power-i-1));	
 	}
     return fval;	
 }
 
 // this function generates the input values
 // in this case, this will be some sawtooth wave
 void gen_data(int size)
 {
 	for(int i=0;i<size;i++)
     {
         data_in_r[i] = (double)i; 
 		data_in_i[i] = 0;
     }	
 	
 	for(int i=size;i<p2_length;i++) // padding with zeros, if necessary
     {
         data_in_r[i] = 0; 
 		data_in_i[i] = 0; 		
     }	
 
 }
  
  // this function gives the next power of two
  int next_2pow(int input)
  {
  	what_power = ceil(log((double)input)/log(2.0));
      return (int)pow(2.0, what_power);
  }
  
  // this function gives the last power of two
  int last_2pow(int input)
  {
      return (int)pow(2.0, floor(log((double)input)/log(2.0)));
  }
  
  // W Twiddle factor - real part
  // 
  double WLn_r(unsigned int L, unsigned int n)  // W is complex, we handle it with a real (_r) and with an imaginary (_i)
  {
  	return ( (double)cos( 2*M_PI*n /L ));
  }
  
  // W Twiddle factor - imaginary part
  double WLn_i(unsigned int L, unsigned int n)
  {
  	return ( (double)sin( 2*M_PI*n /L ));
  }
 
//These functions are not affected by parallelizing the code, except one function which we will use for determining the number of workers. This function is the last_2pow() function and it calculates the lower neighbouring power of two of the input. We will get back to this later on.

//The next change is in the first rows of the main() function.



 int main()
 {
int  namelen;  
  
 	double s_time; // MPI timing
 	// Initialize the MPI environment - we can use the MPI functions and world until MPI_Finalize
     MPI_Init(NULL, NULL);
 	
 	s_time=MPI_Wtime();
 	// number of processes
     MPI_Comm_size(MPI_COMM_WORLD, &world_size);
     // Get the rank of the process
     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
     char processor_name[MPI_MAX_PROCESSOR_NAME];
     MPI_Get_processor_name(processor_name,&namelen); 

     //fprintf(stderr,"Process %d on %s\n", myid, processor_name); 
     cout<<"Process"<<world_rank<<" on "<<processor_name;
 	// we set the required input size to power of two
 	p2_length=next_2pow(p_length); // next power 2 length 
 	bitlength = log(p2_length) / log(2); // bit length  - for monitoring purposes
//We define a double variable named s_time and use it for MPI timing with the MPI_Wtime() function. This later gives a time in seconds elapsed since an arbitrary time in the past. For us, differences between such values will be of use. After initializing the MPI environment (MPI_Init()), we enumerate the workers in our default MPI communicator with MPI_Comm_size() and put the number into world_size. The next step is to identify ourselves and get our rank with the MPI_Comm_rank() function and store it to world_rank. The next two lines are familiar from the serial version, we set the input number to power of two and set some monitoring parameters.

 
	// we have to split the work among workers
	// we wan't power of 2 workers
	int ws_temp = world_size;
	// should we have less data than workers..
	if ( p2_length < world_size )
	{
		ws_temp = p2_length;
		// this will be always power of two, since our p2_length is such
	}
	else
	{
		ws_temp = last_2pow(world_size);  // we use power of 2 number of workers
	}
 
    // unnecessary workers check
	if(world_rank > ws_temp-1)
	{
		// we are not needed, so we quit
		cout << world_rank << ": ---- Bye! ---- " << endl;
		MPI_Finalize();
		exit(1);
	}
//At this point, we have to prepare the data distribution. In this example, we will hand out equal portions of data to every working node. If the amount of data is less than the number of workers, we select data number of workers active, thus we set the variable ws_temp to p2_length. If we have more data than workers, which is mostly the case, we limit the number of workers to the maximal power of two possible. For example, if we have 5 workers and 32 input data, we use only 4 workers. Why do we want to use power of two number of workers? Because we can't do anything else. If we have power of two number of input data and we want to divide it to equal parts, the number of parts will be power of two and the number of data in one part will be a power of two as well. Since we don't need the unnecessary workers, we just let them say "Bye", finalize the MPI environment and exit the program. Although this move works in our example, we have to pay attention to one issue regarding this. If some of the initial members exit, we will loose all MPI functionality regarding collective communication. A workaround for that will be discussed at the end of this subsection.

//The next block is similar to the serial version, except some data range calculations at the beginning.

 	
 
 	 // we divide the number of input data with the number of workers 
 	unsigned int dataslice = p2_length / ws_temp; 
 	FROM = dataslice*world_rank;
 	TO = FROM + dataslice;
 	
   // information about calculated values
   if (world_rank == 0)
   {
     cout << "p_length " << p_length << endl;
     cout << "p2_length " << p2_length << endl;
     cout << "bitlength " << bitlength << endl;
     cout << "what_power " << what_power << endl;
     cout << "No. of workers: " << ws_temp << endl;	
   }
   
  // allocate memory for input data REAL
   data_in_r = new double[p2_length];
 
   // allocate memory for input data IMAGINARY
   data_in_i = new double[p2_length];
   
   // allocate memory for output data (real)
   data_out_r = new double[p2_length];
   
   // allocate memory for output data (im)
   data_out_i = new double[p2_length];
   
 // allocate memory for buffer data (working buffer)
   buffer_r = new double[p2_length];
   
   buffer_i = new double[p2_length];
   
 // allocate memory for reverse binary sorted data (working buffer)
   rb_buffer_r = new double[p2_length];
   
   rb_buffer_i = new double[p2_length];
   
   // generate input data as a saw form
   gen_data(p_length);  // we generate data of p_length length, if it's not a power of two, it will be padded with zeros in the function
   // put data in index-bit reverse order to data_out 
   
   
   /* only a short debug to look at the reverse 
180:   for(int i=0;i<p2_length;i++)
181:   {
182:     std::bitset<32> x(i);
183: 	unsigned int rev = revbit(i);
184: 	std::bitset<32> y(rev);
185: 	// cout << "szam: " << i << " (" << x << ") reverse: " << rev << " (" << y << ") \r\n";
186:   }
187:   */
  // we will start the loops here.
  // we work with complex numbers, so we have a real (_r) and an imaginary (_i) part as well
   
   unsigned int puf_l = p2_length;  // working variable - starting length of input data - it will be changed during the algorithm (halved and halved..) 
   unsigned int divider = 0;  // divider variable, it will be doubled as we step forward with the series of Even Odd recombination levels
   unsigned int ind=0 ;	 // temp variable for copying data
   
   // at first, we generate a reverse binary ordered working buffer from the input. It's worth to do it once and just copy it into the working buffer whenever it's necessary
  // copy input data into buffer in reverse binary order
 // input data is complex, hence the real and imaginary part
  for(unsigned int i=0;i<p2_length;i++)
   {
 		ind = revbit(i);
 		rb_buffer_r[i]=data_in_r[ind];
 		rb_buffer_i[i]=data_in_i[ind];
   }
//First, we set the amount of data for one worker in dataslice and every worker determines its range of work (FROM, TO). In the next few lines, we define some variables, allocate memory for our lists, generate the input data and the reverse binary order data series just as in the serial version.


  
  // the FROM and TO variables are set for every worker different
  for(unsigned int n=FROM;n<TO;n++)   // F(n)
  {
		// copy reverse binary pre-ordered buffer into working buffer, we have to start from this state at every "n"-step
		for(unsigned int i=0;i<p2_length;i++)
		{
				buffer_r[i]=rb_buffer_r[i];
				buffer_i[i]=rb_buffer_i[i];
		}    
		divider=1;				// we will use this variable to indicate "N"
		puf_l = p2_length;    // length starts always from the input length and will change inside the next loop
		for(unsigned int b=0;b<what_power;b++)  // b is running from 0 to logN levels
		{
		    divider = divider*2;  // as we step forward, we have to double the value of N
			for(unsigned int k=0;k<puf_l;k++) 
			{
   			   buffer_r[k] = buffer_r[2*k] + buffer_r[2*k+1]*WLn_r(divider ,n)  - buffer_i[2*k+1]*WLn_i(divider, n);   
			   // a+bi    x      c+di   = (ac - bd) + i(ad + bc)		 
   			   buffer_i[k] = buffer_i[2*k]  + buffer_r[2*k+1]*WLn_i(divider ,n) + buffer_i[2*k+1]*WLn_r(divider, n);    
			   // m+ni     +   (a+bi x c+di)   even + odd * w   -->   r:   m+(ac - bd)      i:  n+(ad + bc)
		    }
			puf_l = puf_l / 2; 
		}
		data_out_r[n]=buffer_r[0];
		data_out_i[n]=buffer_i[0];
		
		// we look at the input and output data
		// cout << "-- Input " << n << ": " << data_in_r[n] << "," << data_in_i[n] << " -- output: " << data_out_r[n] << "," << data_out_i[n] << endl;
	}
//The outer loop, which processes the input data is running on a definite interval for every worker. The inside of this loop doesn't change compared to the serial code. Since every worker has the whole dataset at the beginning and the algorithm doesn't rely on previous results, we can run the algorithm on every worker parallel. After doing the FFT on the given input range, every worker has its own result chunk. This has to be collected to get the whole result at one process.


 	
    // we have to collect all data
    if(ws_temp>1)  // we collect only, if there are more than one processes
 	{   
 		for(int p=1;p<ws_temp;p++) // we have to collect ws_temp-1 slices
 		{
 			if(world_rank == p) // if it's my turn, as sender, I send.
 			{
 			    // cout << "I'm sending (me: " << world_rank << ") data from " << FROM << " to " << TO << ". dataslice: " << dataslice << endl;
 				MPI_Send(&data_out_r[FROM], dataslice, MPI_DOUBLE, 0 , 900 , MPI_COMM_WORLD); 
 				MPI_Send(&data_out_i[FROM], dataslice, MPI_DOUBLE, 0 , 900 , MPI_COMM_WORLD); 
 				
 			}
 			else if ( world_rank == 0) // if we are process 0, we collect the sended data
			{
				   // what shall rFROM and rTO be? What lines do we have to collect. We have to calculate it from the remote process' rank (p) 
					rFROM = dataslice*p;
					rTO = dataslice*(p+1);
							
				// cout << "I receive (me: " << world_rank << ") the lines from: " << p << ". Interval: " << rFROM << " - " << rTO << endl;
				MPI_Recv(&data_out_r[rFROM], dataslice, MPI_DOUBLE, p , 900 , MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
				MPI_Recv(&data_out_i[rFROM], dataslice, MPI_DOUBLE, p , 900 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
			} 		}		 	
		}
 	if(world_rank == 0)
 	{
 		cout << world_rank << " collecting done! - MPI_time: " << MPI_Wtime()-s_time << endl;
 	}
//First, we check whether more than one workers are running. If yes, we have to collect data from ws_temp-1 co-workers. It is always a good idea to collect data at the first worker with rank zero. This rank is always present. If there are more workers, then every worker with a greater rank than zero has to send (MPI_Send()) the real (data_out_r[]) and imaginary (data_out_i[]) part of its data chunk to process 0. Parallel to this, the zero rank worker has to receive (MPI_Recv()) these data chunks (dataslice) into the right place of the real and imaginary arrays. Although the right place has to be calculated before. Since we know the rank of the other workers we are waiting data from, we can easily calculate the starting (rFROM) and ending (rTO) points of the data chunks we are getting from the remote processes.

	// some Monitoring output to look what data did we calculate and collect
	if(world_rank == 0)
	{ 	  for(int n=0;n<p2_length;n++) 	  cout << "-- Input " << n << ": " << data_in_r[n] + data_in_i[n] << " -- output: " << data_out_r[n] - data_out_i[n] << endl; 	} 
   // we loose data after this lines, any output should be done here
   delete[] data_out_i; // free up memory
   delete[] data_out_r;
   delete[] data_in_i;
   delete[] data_in_r;
   delete[] buffer_r;
   delete[] buffer_i;
   delete[] rb_buffer_r;
   delete[] rb_buffer_i;
   
 	// MPI end - we loose contact with the other nodes at this point
 	MPI_Finalize();
   
   return 0;
 }