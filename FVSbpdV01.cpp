/* 
********************************************************************************
FVSbpdV01.cpp

DESCRIPTION:
Belief-propagation guided decimation (BPD) as a solver for the maximum feedback
vertex set problem. This program is applicable on a single graph instance. The
imput graph file has the following form:

N   M                % first row specifies N (vertex number) and M (edge number)
i_1  j_1                                            % undirected edge (i_1, j_1)
i_2  j_2                                            % undirected edge (i_2, j_2)
.    .
.    .
.    .
i_M  j_M                                            % undirected edge (i_M, j_M)

The program reads only the first EdgeNumber edges from the imput graph, with
EdgeNumber being explicitly specified (EdgeNumber <= M, of course).

Each vertex has two states (unoccupied, 0; occupied, 1). An unoccupied vertex
belongs to the constructed feedback vertex set S, while an occupied vertex does
not belong to S.

Initially all the vertices are occupied and active. At each round of decimation,
some of the remaining active vertices are fixed to be un-occupied sequentially.
After a vertex is fixed to be un-occupied (and it becomes inactive), the
remaining active subgraph is simplified (during which more vertices become
inactive). During the simplication process, all the leaf vertices in the active
subgraph are removed recursively (and they are fixed to be occupied).

For details about the replica-symmetric mean field theory and the spin glass
model on which this algorithm was based, please consult the following
references:
[1] Hai-Jun Zhou, "Spin Glasses and Message Passing" (Science Press, Beijing,
2015), chapter 7, pages 218--238.
[2] Hai-Jun Zhou, "Spin glass approach to the feedback vertex set problem",
European Physical Journal B 86: 455 (2013).

To generate the executive file, you can simply compile as
* c++ -O3 -o bpdfvs.exe FVSbpdV01.cpp

!!!please make sure some of the key parameters, such as input graph name, number
of edges, number of vertices, output file names, as appropriately specified in 
the FVSbpdV01.cpp file.

After you successfully compiled the code, you can then run the algorithm as
* bpdfvs.exe

Good luck!

LOG:
06.09.2015--08.09.2015: FVSbpdV01.cpp revision. Make the algorithm works even
if some of the vertices have extremely many attached edges.
06.09.2015: FVSbpdV01.cpp (copied from FVSbpdv00.cpp)
29.06.2013--05.07.2013: FVSbpdv00.cpp

PROGRAMMER:
Hai-Jun Zhou
Institute of Theoretical Physics, Chinese Academy of Sciences
Zhong-Guan-Cun East Road 55, Beijing 100190
email: zhouhj@itp.ac.cn
webpage: power.itp.ac.cn/~zhouhj/
********************************************************************************
*/

#include <exception>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <valarray>

//#include "/Users/zhouhj/Programs/zhjrandom.h"      //a random number generator
#include "zhjrandom.h"                             //a random number generator

using namespace std;

/*---               random real uniformly distributed in [0,1)            ---*/
double u01prn(ZHJRANDOMv3 *rd)
{
  return rd->rdflt();
}

struct OrderedIntPair
{
  int first;
  int second;
  OrderedIntPair(void);
  OrderedIntPair(int, int);
};

OrderedIntPair::OrderedIntPair(void)
{
  first=0;
  second=0;
  return;
}

OrderedIntPair::OrderedIntPair(int a, int b)
{
  if(a<=b)
    {
      first  = a;
      second = b;
    }
  else
    {
      first  = b;
      second = a;
    }
  return ;
}

bool operator<(OrderedIntPair a, OrderedIntPair b)
{
  if(a.first < b.first)
    return true;
  else if(a.first > b.first)
    return false;
  else                                                       //a.first = b.first
    {
      if(a.second < b.second)
	return true;
      else
	return false;
    }
}

bool operator==(OrderedIntPair a, OrderedIntPair b)
{
  return (a.first == b.first) && (a.second == b.second);
}

struct message                                    //cavity message from a vertex
{
  struct vstruct *v_ptr;                         //pointer to the sending vertex
  double q_0;                                  //probability of taking state A=0
  double q_root;                                 //probability of being the root
  message(void);
  message(double, double);
};

message::message(void)
{
  v_ptr  = 0;
  q_0    = 0;
  q_root = 0;
  return;
}

message::message(double a, double b)
{
  v_ptr  = 0;
  q_0    = a;
  q_root = b;
  return;
}

struct mpointer                                             //pointer to message
{
  struct message *m_ptr;
  mpointer(void);
};

mpointer::mpointer(void)
{
  m_ptr=0;
  return;
}

struct vstruct                                                   //vertex struct
{
  int index;                                 //index of vertex, positive integer
  int degree;                                   //number of neighboring vertices
  int active_degree;                     //number of active neighboring vertices
  bool occupied;           //vertex state A (there are degree+2 possible values)
  bool active;  //=true (not being removed from graph); =true (has been deleted)
  double q_0;                                                //empty probability

  struct message *im_ptr;                      //start position of input message
  struct mpointer *omptr_ptr;    //start position of output message address list
  vstruct(void);
};

vstruct::vstruct(void)
{
  index=0;
  degree=0;
  active_degree=0;
  occupied=true;                           //initially all vertices are occupied
  active=true;                               //initially all vertices are active
  q_0=0;
  im_ptr=0;
  omptr_ptr=0;
  return;
}

class FVS                                                  //feedback vertex set
{
public:
  FVS(ZHJRANDOMv3* );                                              //constructor
  void SetX(double);                              //set re-weighting parameter X
  void SetDampingFactor(double);                            //set damping factor
  bool Graph(const string&, int);                //read graph connection pattern
  
  void Initialization(void);                             //initialize population
  bool ReadPopulation(ifstream& );             //initialize population from file
  bool CheckFVS(const string&);      //check whether occupation pattern is a FVS
  void ReportPopulation(ofstream& );        //record population for future usage
  
  bool BeliefPropagation(double, int);                     //population updating
  void EmptyProbability(void);        //empty probability for each active vertex
  void Thermodynamics(const string& );                //thermodynamic quantities
  int  Fix0(void);           //fix some variables to be un-occupied and simplify
  void Simplify(struct vstruct* );         //simplify the graph by removing leaf
  
private:
  int VertexNumber;                                   //total number of vertices
  int ActiveVertexNumber;                      //total number of active vertices
  int ActiveVertexNumber0;    //total # of active vertices before this BPD round
  
  int EdgeNumber;                                        //total number of edges
  int MaxDegree;                            //maximal vertex degree in the graph

  int MaxFixNumber;              //maximal number of fixed vertices in one round
  int MinRange;                  //minimal range of empty_prob of fixed vertices
  float FixFraction;            /* fraction of active vertices to be fixed to be
				   un-occupied in each round of decimation */
  double X;                                             //re-weighting parameter
  double Weight0;                                                     //=exp(-X)
  double DampingFactor;                     //damping factor of message updating
  
  ZHJRANDOMv3 *PRNG;                                   //random number generator
  
  valarray<vstruct> Vertex;
  valarray<message> InMessage;
  valarray<mpointer> OutMessageAddress;
 
  valarray<int> CandidateVertices;      //list of candidate vertices to be fixed
  valarray<int> CandidateSize;      //num of candidates in each empty_prob range
  valarray<int> Permutation;          //array used in random sequential updating
  valarray<double> WeightA;        //auxiliary array needed for message updating
  valarray<double> WeightB;        //auxiliary array needed for message updating
  valarray<double> WeightC;        //auxiliary array needed for message updating
  valarray<bool>   IsActive;       //auxiliary array needed for message updating
 
  void UpdateMessage(struct vstruct *, double&);
};

int main(int argc, char ** argv)
{
  //random number generator initialization
  int rdseed=54716791;                  //you can set this seed to another value
  ZHJRANDOMv3 rdgenerator(rdseed);
  int prerun=10000000;                         //you can set it to another value
  for(int i=0; i<prerun; ++i)
    rdgenerator.rdflt();
  
  FVS system( &rdgenerator);
  
  //  read graph, please check these parameters according to your graph instance
  int VertexNumber = 100000;                   //number of vertices in the graph
  int EdgeNumber   = 500000;    //number of edges you want to read from the file
  string graphname="ERn100kM1m.g100";                        //graph file's name

  bool succeed=system.Graph(graphname, EdgeNumber);
  if(succeed==false)
    return -1;
  
  //initialize cavity messages in a most random way
  system.Initialization();
  
  /*
    ifstream inputf("temp.dat");
    system.ReadPopulation(inputf);                  // cavity messages from file
  */
  
  //set BPD parameters
  double X=12.0e0;
  system.SetX(X);
  double DampingFactor=0.9e0;
  system.SetDampingFactor(DampingFactor);
  
  //initial BP iterations
  float epsilon = 1.0e-15;
  int eq_iterations = 500;     //please don't make it too small (i.e., >= 100)
  system.BeliefPropagation(epsilon, eq_iterations);
  
  string thermfile="ERn100kM1mg100FreeEnergy.dat";
  system.Thermodynamics(thermfile);                      //compute free energy
  
  /*
    ofstream outputf("temp.dat");
    system.ReportPopulation(outputf);           //output cavity messages to file
  */
  
  //BPD
  epsilon = 1.0e-7;
  int iterations = 50;  /*the final result not sensitive to this parameter, as
			  long as it is not too small, e.g., iterations>10 */
  int ActiveVertexNumber=0;
  do
    {
      system.BeliefPropagation(epsilon, iterations);
      system.EmptyProbability();
      
      ActiveVertexNumber=system.Fix0();
    }
  while(ActiveVertexNumber>0);
  
  //            report feedback vertex set, please change the name as you prefer
  string FVSfile="ERn100kM1mg100.FVS";
  if( system.CheckFVS(FVSfile) == true)
    {
      cout<<"FVS constructed.\n";
      return 1;
    }
  else
    {
      cerr<<"Not a feedback vertex set.\n";
      return -1;
    }
}

/*---                       constructor of FVS cluster                     ---*/
FVS::FVS(ZHJRANDOMv3* rd)
{
  PRNG = rd;                                           //random number generator
  FixFraction = 0.01;                     //fix 1 percent of the active vertices
  return;
}

/* -                            Read graph from file                           -
   gname: the input graph name.
   enumber: read the first enumber edges only.
*/
bool FVS::Graph(const string& gname, int enumber)
{
  ifstream graphf(gname.c_str());
  if( !graphf.good() )
    {
      cerr<<"Graph probably non-existant.\n";
      return false;
    }
  
  //first read of input graph
  graphf >>VertexNumber
	 >>EdgeNumber;
  
  if(EdgeNumber<enumber)
    {
      cerr<<"No so many edges in the graph.\n";
      graphf.close();
      return false;
    }
  
  EdgeNumber=enumber;                  //only the first enum edges are read into
  
  try { Vertex.resize(VertexNumber+1); } catch(bad_alloc)
    {
      cerr<<"Vertex construction failed.\n";
      graphf.close();
      return false;
    }
  
  bool succeed=true;
  set<OrderedIntPair> EdgeSet;
  for(int eindex=0; eindex<EdgeNumber && succeed; ++eindex)
    {
      int v1, v2;
      graphf >>v1
	     >>v2;
      if(v1==v2 || v1==0 || v1 >VertexNumber || v2==0 || v2 >VertexNumber)
	{
	  cerr<<"Graph incorrect at line "<<eindex+1<<endl;
	  succeed=false;
	}
      else if(EdgeSet.find(OrderedIntPair(v1,v2) ) != EdgeSet.end() )
	{
	  cerr<<"Multiple edges.\n";
	  succeed=false;
	}
      else
	{
	  EdgeSet.insert(OrderedIntPair(v1,v2));
	  ++(Vertex[v1].degree);
	  ++(Vertex[v2].degree);
	}
    }
  graphf.close();
  if(succeed==false)
    return false;
 
  EdgeSet.clear();
  try { InMessage.resize( 2*EdgeNumber ); } catch(bad_alloc)
    { 
      cerr<<"InMessage construction failed.\n";
      return false;
    }
  try { OutMessageAddress.resize( 2*EdgeNumber ); } catch(bad_alloc)
    {
      cerr<<"OutMessageAddress construction failed.\n";
      return false;
    }
  
  int position=0;
  MaxDegree=0;
  struct vstruct *v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v)
    {
      v_ptr->index =v;
      v_ptr->occupied=true;
      v_ptr->active =true;
      v_ptr->active_degree = v_ptr->degree;
      
      v_ptr->im_ptr = &InMessage[position];
      v_ptr->omptr_ptr = &OutMessageAddress[position];
      position += v_ptr->degree;
      if(v_ptr->degree>MaxDegree)
	MaxDegree = v_ptr->degree;
      v_ptr->degree=0;
      ++v_ptr;
    }
  
  //second read of input graph
  graphf.open(gname.c_str());
  graphf >>VertexNumber
	 >>EdgeNumber;
  EdgeNumber=enumber;
  
  struct mpointer *omptr_ptr=&OutMessageAddress[0];
  for(int eindex=0; eindex<EdgeNumber; ++eindex)
    {
      int v1,v2;
      graphf >>v1
	     >>v2;
      omptr_ptr=Vertex[v1].omptr_ptr + Vertex[v1].degree;
      omptr_ptr->m_ptr = Vertex[v2].im_ptr + Vertex[v2].degree;
      omptr_ptr->m_ptr->v_ptr = &Vertex[v1];
      
      omptr_ptr=Vertex[v2].omptr_ptr + Vertex[v2].degree;
      omptr_ptr->m_ptr = Vertex[v1].im_ptr + Vertex[v1].degree;
      omptr_ptr->m_ptr->v_ptr = &Vertex[v2];
      
      ++(Vertex[v1].degree);
      ++(Vertex[v2].degree);
    }
  graphf.close();
  
  cout << "Graph: N= "<<VertexNumber <<",  M= "
       << EdgeNumber   <<",  Mean degree= "
       << (2.0e0*EdgeNumber)/(1.0e0*VertexNumber)<<endl<<endl;
  
  try { WeightA.resize(MaxDegree); } catch(bad_alloc)
    {
      cerr<<"WeightA construction failed.\n";
      return false;
    }
  try { WeightB.resize(MaxDegree); } catch(bad_alloc)
    {
      cerr<<"WeightB construction failed.\n";
      return false;
    }
  try { WeightC.resize(MaxDegree); } catch(bad_alloc)
    {
      cerr<<"WeightC construction failed.\n";
      return false;
    }
  try { IsActive.resize(MaxDegree); } catch(bad_alloc)
    {
      cerr<<"IsActive construction failed.\n";
      return false;
    }

  ActiveVertexNumber=VertexNumber;
  v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v)
    {
      if(v_ptr->active)
	{
	  if(v_ptr->active_degree <= 1 )
	    {
	      v_ptr->active = false;
	      v_ptr->occupied = true;                          //being occupied
	      --ActiveVertexNumber;
	      Simplify(v_ptr);
	    }
	}
      ++v_ptr;
    }

  if(ActiveVertexNumber==0)
    {
      cerr<<"The graph has no loops. Done! \n";
      return false;
    }
  
  try { Permutation.resize(ActiveVertexNumber); } catch(bad_alloc) 
    {
      cerr<<"Permutation construction failed.\n";
      return false;
    }
  v_ptr=&Vertex[1];
  ActiveVertexNumber=0;
  for(int v=1; v<=VertexNumber; ++v)
    {
      if(v_ptr->active == true)
	{
	  Permutation[ActiveVertexNumber]=v_ptr->index;
	  ++ActiveVertexNumber;
	}
      ++v_ptr;
    }
  ActiveVertexNumber0=ActiveVertexNumber;

  MaxFixNumber=static_cast<int>(ActiveVertexNumber*FixFraction);
  if(MaxFixNumber==0)
    MaxFixNumber = 1;
  
  //the empty probability q_0 is divided into bins of width 0.01
  try { CandidateVertices.resize( MaxFixNumber * 101 ); } catch( bad_alloc )
    {
      cerr<<"CandidateVertices construction failed.\n";
      return false;
    }
  CandidateSize.resize(101); 
  
  MinRange=0;
  
  return true;
}

/*---                set the value of X                                   ---*/
void FVS::SetX(double xval)
{
  X = xval;                                             //re-weighting parameter
  Weight0=exp(-X);
  return;
}

/*---                set the value of DampingFactor                       ---*/
void FVS::SetDampingFactor(double dp)
{
  DampingFactor = dp;                                           //damping factor
  if(DampingFactor >= 1.0e0)
    DampingFactor=1.0e0;
  else if(DampingFactor<=0.0001e0)
    DampingFactor=0.0001e0;
  return;
}

/*---                      initialize population                          ---*/
void FVS::Initialization(void)
{
  struct mpointer *omptr_ptr = &OutMessageAddress[0];
  for(int v=1; v<=VertexNumber; ++v)
    {
      int degree=Vertex[v].degree;
      for(int d=0; d<degree; ++d)
	{
	  omptr_ptr->m_ptr->q_0 = 1.0e0/(1.0e0*degree+2.0e0);
	  omptr_ptr->m_ptr->q_root = 1.0e0/(1.0e0*degree+2.0e0);
	  ++omptr_ptr;
	}                                                 //uniform distribution
    }
  return;
}

/*---                initialize population from file                      ---*/
bool FVS::ReadPopulation(ifstream& ipopfile)
{
  int enumber;
  
  ipopfile>>enumber;
  if(enumber != EdgeNumber)
    {
      cerr<<"Edge number does not match.\n";
      return false;
    }
  
  double q_0, q_root;
  struct message *im_ptr = &InMessage[0];
  for(int index=0; index<2*EdgeNumber; ++index)
    {
      ipopfile >> q_0 >> q_root;
      im_ptr->q_0 = q_0;
      im_ptr->q_root = q_root;
      ++im_ptr;
    }
  return true;
}

/*---                output population to file                            ---*/
void FVS::ReportPopulation(ofstream& opopfile)
{
  opopfile << EdgeNumber <<endl<<endl;
  
  struct message *im_ptr = &InMessage[0];
  int degree;
  for(int v=1; v<=VertexNumber; ++v)
    {
      im_ptr = Vertex[v].im_ptr;
      degree = Vertex[v].degree;
      for(int d=0; d<degree; ++d)
	{
	  opopfile << im_ptr->q_0 <<'\t'
		   << im_ptr->q_root << endl;
	  ++im_ptr;
	}
      opopfile<<endl;
    }
  
  return;
}

/*                 simplify after fixing variable                             */
void FVS::Simplify(struct vstruct *v_ptr)
{
  struct message *im_ptr=v_ptr->im_ptr;
  for(int d=0; d < v_ptr->degree; ++d)
    //                               if central vertex has at least one neighbor
    {
      if( --(im_ptr->v_ptr->active_degree) <= 1)
	//                            vertex now has zero or one active neighbor
	{
	  if(im_ptr->v_ptr->active)
	    //                                       and vertex is itself active
	    {
	      im_ptr->v_ptr->active   = false;
	      im_ptr->v_ptr->occupied = true;
	      --ActiveVertexNumber;
	      Simplify(im_ptr->v_ptr);
	    }
	}
      ++im_ptr;
    }
  return;
}

/* -                 Belief propagation                                     - */
bool FVS::BeliefPropagation(double error, int count)
{
  int iter=0;
  double max_error=0;
  
  for(int quant=ActiveVertexNumber0-1; quant>=0; --quant)
    {
      if(Vertex[ Permutation[quant] ].active == false)
	{
	  --ActiveVertexNumber0;
	  Permutation[quant]=Permutation[ActiveVertexNumber0];
	}
    }
  
  if(ActiveVertexNumber0 != ActiveVertexNumber) 
    {
      cerr<<"Unexpected mistake.\n";
      return false;
    }
  
  do
    {
      max_error=0;
      for(int quant=ActiveVertexNumber; quant>0; --quant)
	{
	  int iii = static_cast<int>(quant * u01prn(PRNG) );
	  int v = Permutation[iii];
	  Permutation[iii] = Permutation[quant-1];
	  Permutation[quant-1] = v;
	  UpdateMessage(&Vertex[v], max_error);
	}
      if((iter % 100) == 0)
	{ 
	  cerr<<' '<<max_error<<' ';
	  cerr.flush();
	}
      ++iter;
    }
  while (max_error>error && iter<count);
  
  if(max_error<=error)
    {
      cerr<<' '<<max_error<<"  :-)\n";
      return true;
    }
  else
    {
      cerr<<' '<<max_error<<"  :-(\n";
      return false;
    }
}

/* -          Update the output messages from a vertex                         -
   It is assumed that central vertex is active and has at least two active
   neighbors, namely v_ptr->active = true,  v_ptr->active_degree >= 2         */
void FVS::UpdateMessage( struct vstruct *v_ptr, double& maxdiff)
{
  int degree=v_ptr->degree;
  
  struct message *im_ptr = v_ptr->im_ptr;
  for(int j=0; j<degree; ++j)
    {
      if(im_ptr->v_ptr->active)
	{
	  IsActive[j]=true;
	  WeightA[j]=1.0e0; /* probability all neighbors of central vertex i 
			       (except vertex j) are empty or are root in
			       cavity graph */
	  WeightB[j]=0.0e0; /* probability that one neighbor (k not equal to j)
			       of central vertex i are not empty in the cavity
			       graph, all the other neighbors (except j) of
			       central vertex i are either empty or root in
			       the cavity graph */
	  WeightC[j]=Weight0;
	}
      else
	IsActive[j]=false;
      
      ++im_ptr;
    }
  
  im_ptr = v_ptr->im_ptr;
  for(int j=0; j<degree; ++j)
    {
      if(IsActive[j])
	{
	  double q_0     = im_ptr->q_0;
	  double q_root  = im_ptr->q_root;
	  double q_0root = q_0+q_root;
	  
	  for(int k=0; k<degree; ++k)
	    {
	      if((k != j) && IsActive[k])
		{
		  WeightB[k]  = WeightA[k]*(1.0e0-q_0) + WeightB[k]*q_0root;
		  WeightA[k] *= q_0root;
		  
		  //                the following lines were added on 07.09.2015
		  double maxval=max(WeightA[k], WeightB[k]);
		  if(maxval<WeightC[k])
		    maxval=WeightC[k];
		  WeightA[k] /= maxval;                        //avoid underflow
		  WeightB[k] /= maxval;                        //avoid underflow
		  WeightC[k] /= maxval;                        //avoid underflow
		}
	    }
	}
      ++im_ptr;
    }
  
  struct mpointer *omptr_ptr = v_ptr->omptr_ptr;
  for(int j=0; j<degree; ++j)
    {
      if(IsActive[j])
	{
	  double norm=WeightA[j]+WeightB[j]+WeightC[j];
	  double diff_q_0    = WeightC[j]/norm - omptr_ptr->m_ptr->q_0 ;
	  double diff_q_root = WeightA[j]/norm - omptr_ptr->m_ptr->q_root;
	  double diff = abs(diff_q_0)+abs(diff_q_root);
	  if(diff>maxdiff)
	    maxdiff=diff;
	  omptr_ptr->m_ptr->q_0    += DampingFactor * diff_q_0;
	  omptr_ptr->m_ptr->q_root += DampingFactor * diff_q_root;
	}
      ++omptr_ptr;
    }
  
  return;
}

/*---                       Empty probability of each vertex              ---*/
void FVS::EmptyProbability(void)
{
  if(ActiveVertexNumber == 0)
    return ;
  
  double q_0, q_root, q_0root, norm, weight_a, weight_b, weight_c;
  
  struct vstruct *v_ptr=&Vertex[1];
  
  MaxFixNumber=static_cast<int>(ActiveVertexNumber*FixFraction);
  if(MaxFixNumber==0)
    MaxFixNumber = 1;
  
  for(int r=0; r<=100; ++r) CandidateSize[r] = 0;   /* number of candidate
						       vertices in each q_0 bin
						       of width 0.01 */
  MinRange=0;   /* the maximal value of r0 such that 
		   sum_{r>=r0} CandidateSize[r] >= MaxFixNumber is satisfied */
  
  /*         thermodynamic quantities of edges                               */
  struct message *im_ptr = &InMessage[0];
  for(int v_index=0; v_index < ActiveVertexNumber; ++v_index)
    {
      v_ptr = &Vertex[ Permutation[v_index] ];
      int degree=v_ptr->degree;
      
      /*  the central vertex is assumed to have at least two or more active
	  neighbors. If it has zero or only one active neighbor, it should have
	  been removed from the active subgraph */
      
      weight_a=1.0e0; /*          probability all neighbors of central vertex i
				  are empty or are root in cavity graph */
      weight_b=0.0e0; /*     probability that one neighbor (j)of central vertex
			     i are not empty in the cavity graph, all the other
			     neighbors of central vertex i are either empty or
			     root in the cavity graph */
      weight_c=Weight0;
      
      im_ptr = v_ptr->im_ptr;
      for(int j=0; j<degree; ++j)
	{
	  if(im_ptr->v_ptr->active)         // the neighboring vertex is active
	    {
	      q_0     = im_ptr->q_0;
	      q_root  = im_ptr->q_root;
	      q_0root = q_0+q_root;
	      weight_b  = weight_a * (1.0e0-q_0)+weight_b * q_0root;
	      weight_a *= q_0root;
	      
	      //          following lines added on 07.09.2015 to avoid underflow
	      double maxval=max(weight_a, weight_b);
	      if(maxval<weight_c)
		maxval=weight_c;
	      weight_a /= maxval;
	      weight_b /= maxval;
	      weight_c /= maxval;
	    }
	  ++im_ptr;
	}
      norm = weight_a+weight_b+weight_c;
      q_0  = weight_c/norm;
      v_ptr->q_0 = q_0;
      
      int rrr = static_cast<int>(q_0*100);
      if(rrr>=MinRange)
	{
	  if(CandidateSize[rrr]<MaxFixNumber)
	    {
	      CandidateVertices[rrr*MaxFixNumber + CandidateSize[rrr] ] 
		= v_ptr->index;
	      ++CandidateSize[rrr];
	    }
	  else
	    MinRange=rrr;
	}
    }
  
  return;
}

/*---                        thermodynamic quantities                      ---*/
void FVS::Thermodynamics(const string& ofile)
{
  if(ActiveVertexNumber == 0)
    return;
  
  int degree;
  double q_0, q_root, q_0out, q_rootout, q_0root, norm, weight_a, weight_b,
    weight_c, logval;
  
  double phi_vtx=0;                        //vertex contribution to free entropy
  double phi_edge=0;                         //edge contribution to free entropy
  double rho_vtx=0;                                  //vertex occupation density
  
  /*         thermodynamic quantities of edges                               */
  struct message *im_ptr = &InMessage[0];
  struct mpointer *omptr_ptr = &OutMessageAddress[0];
  
  for(int v_index=0; v_index < ActiveVertexNumber; ++v_index)
    {
      struct vstruct *v_ptr= &Vertex[ Permutation[v_index] ];
      int degree=v_ptr->degree;
      
      /*  the central vertex is assumed to have at least two or more active
	  neighbors. If it has zero or only one active neighbor, it should have
	  been removed from the active subgraph */
      
      weight_a=1.0e0; /*          probability all neighbors of central vertex i
				  are empty or are root in cavity graph */
      weight_b=0.0e0; /*     probability that one neighbor (j)of central vertex
			     i are not empty in the cavity graph, all the other
			     neighbors of central vertex i are either empty or
			     root in the cavity graph */
      weight_c=Weight0;
      
      logval=0;
      im_ptr = v_ptr->im_ptr;
      omptr_ptr = v_ptr->omptr_ptr;
      for(int j=0; j<degree; ++j)
	{
	  q_0out = omptr_ptr->m_ptr->q_0;
	  q_rootout=omptr_ptr->m_ptr->q_root;
	  
	  if(im_ptr->v_ptr->active)          // the neighboring vertex is active
	    {
	      q_0     = im_ptr->q_0;
	      q_root  = im_ptr->q_root;
	      q_0root = q_0+q_root;
	      weight_b  = weight_a * (1.0e0-q_0)+weight_b * q_0root;
	      weight_a *= q_0root;
	      
	      //  the following lines added on 07.09.2015 for avoiding underflow
	      double maxval = max(weight_a, weight_b);
	      if(maxval<weight_c)
		maxval=weight_c;
	      weight_a /= maxval;
	      weight_b /= maxval;
	      weight_c /= maxval;
	      logval += log(maxval);

	      phi_edge += log(q_0 + q_0out * (1.0e0-q_0) + q_root * 
			      (1.0e0 - q_0out) + q_rootout * (1.0e0 - q_0) );
	    }
	  ++im_ptr;
	  ++omptr_ptr;
	}
      norm = weight_a+weight_b+weight_c;
      q_0  = weight_c/norm;
      v_ptr->q_0 = q_0;
      
      rho_vtx += 1.0e0-q_0;
      //the next line was modified on 07.09.2015: adding the last term "+logval"
      phi_vtx += X+log(norm)+logval;
    }
  
  phi_vtx  /= ActiveVertexNumber;
  phi_edge /= 2.0e0*ActiveVertexNumber;
  rho_vtx  /= ActiveVertexNumber;
  
  double phi = (phi_vtx - phi_edge)/X;
  ofstream output(ofile.c_str() );
  output<< X << '\t'
	<< ActiveVertexNumber << '\t'
	<< rho_vtx <<'\t'
	<< phi <<'\t'
	<< X * (phi-rho_vtx) << endl;

  cout<<"Estimated FVS relative size at X="<<X<<" is "<<1-rho_vtx<<endl;
  
  return ;
}

/* -- externally fixing some variables to be empty and simplify the system --*/
int FVS::Fix0(void)
{
  int rank = 100; /* q_0 is distributed to 101 bins [0,0.01), [0.01,0.02), ...,
		     [0.99,1), [1,1] */
  
  int num_fixed_empty = 0;  //number of vertices externally fixed to unoccupied
  double mean_emptyprob = 0;        // ... and mean q_0 value of these vertices
  
  int num_examined_vertices = 0;        //number of examined candidate vertices
  struct vstruct *v_ptr;
  
  while(num_examined_vertices < MaxFixNumber && ActiveVertexNumber > 0)
    {
      int size=CandidateSize[rank];
      if(size>0)
	{
	  int *i_ptr = &CandidateVertices[rank * MaxFixNumber];
	  if(num_examined_vertices + size <= MaxFixNumber)
	    {
	      num_examined_vertices += size;
	      for(int s=0; s<size; ++s)
		{
		  v_ptr = &Vertex[ *i_ptr ];
		  if(v_ptr->active)
		    {
		      ++num_fixed_empty;
		      v_ptr->active = false;
		      v_ptr->occupied = false;
		      mean_emptyprob += v_ptr->q_0;
		      --ActiveVertexNumber;
		      
		      Simplify(v_ptr) ;
		    }
		  ++i_ptr;
		}
	    }
	  else
	    {
	      while(num_examined_vertices < MaxFixNumber )
		{
		  ++num_examined_vertices;
		  v_ptr = &Vertex[ *i_ptr ];
		  if(v_ptr->active)
		    {
		      v_ptr->active = false;
		      v_ptr->occupied = false;
		      ++num_fixed_empty;
		      mean_emptyprob += v_ptr->q_0;
		      --ActiveVertexNumber;
		      
		      Simplify(v_ptr) ;
		    }
		  ++i_ptr;
		}
	    }
	}
      --rank;
    }
  
  cout<<" - number of active vertices: "<<ActiveVertexNumber<<endl;

  return ActiveVertexNumber;
}

/*  --- check whether the final occupation pattern corresponds to a FVS   --- */
bool FVS::CheckFVS( const string& filename)
{
  int num_occup =0;
  int num_empty =0;                                                   //FVS size
  int num_active_edge =0;
  
  /* Consider the subgraph induced by all the occupied vertices and the edges
     between pairs of such vertices. 
     - For each vertex in this subgraph, first determine its number of neighbors
     in this subgraph.
     - Then check whether all the edges of this subnetwork can be completely
     deleted by removing vertices of degree one iteratively.  */
  queue<int> LeafVertices;
  struct vstruct *v_ptr = &Vertex[1];
  struct message *im_ptr = &InMessage[0];
  for(int v=1; v<=VertexNumber; ++v)
    {
      if(v_ptr->occupied == false)
	++num_empty;                                     //vertex belongs to FVS
      else
	{
	  ++num_occup;
	  im_ptr = v_ptr->im_ptr;
	  v_ptr->active_degree=0;                   //number of active neighbors
	  for(int d=0; d < v_ptr->degree; ++d)
	    {
	      if(im_ptr->v_ptr->occupied)
		{
		  ++(v_ptr->active_degree);
		  ++num_active_edge;
		}
	      ++im_ptr;
	    }
	  
	  if( v_ptr->active_degree == 1)
	    LeafVertices.push(v_ptr->index);
	}
      ++v_ptr;
    }
  num_active_edge /= 2;
  
  while( !LeafVertices.empty() )
    {
      v_ptr = &Vertex[ LeafVertices.front() ];
      LeafVertices.pop();
      
      im_ptr = v_ptr->im_ptr;
      for(int d=0; d < v_ptr->degree; ++d)
	{
	  if(im_ptr->v_ptr->occupied)
	    {
	      if( --(im_ptr->v_ptr->active_degree) == 1)
		LeafVertices.push(im_ptr->v_ptr->index);
	    }
	  ++im_ptr;
	}
    }
  
  v_ptr = &Vertex[1];
  int TwoCoreSize=0;
  for(int v=1; v<=VertexNumber; ++v)
    {
      if(v_ptr->occupied && v_ptr->active_degree != 0)
	++TwoCoreSize;
      ++v_ptr;
    }
  
  if(TwoCoreSize==0)
    {
      ofstream pfile(filename.c_str() );
      
      pfile <<"feedback-vertex-set-size: "<<num_empty<<endl<<endl;

      cout<< "FVS size = "<<num_empty<<endl;
      
      v_ptr=&Vertex[1];
      for(int v=1; v<=VertexNumber; ++v)
	{
	  if(v_ptr->occupied==false)
	    pfile << v_ptr->index<<endl;
	  ++v_ptr;
	}
      pfile.close();
      
      return true;
    }
  else
    {
      cerr<<"Not a proper FVS. The final two-core size is "<<TwoCoreSize<<endl;
      return false;
    }
}
