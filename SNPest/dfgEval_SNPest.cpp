#include <boost/program_options.hpp>
#include "phy/DfgIO.h"

namespace po = boost::program_options;
using namespace phy;

// SL: I added this function
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

string itoa(int a){
  string res="";
  int remainder;
  if(a == 0) res="0";
  while(a>0){
    remainder=a%10;
    a=int(a/10);
    res=char(remainder+48)+res;
  }
  return res;
}

void takeMinusLog(xvector_t & v, string const & id)
{
  for (unsigned i = 0; i < v.size(); i++) {
    if (v[i] < 0)
      errorAbort("From takeMinusLog: Trying to logarithm of negative number resulting from data line with id: '" + id + "'.");
    v[i] = - log(v[i]);
  }
}


// useful with post probs close to one, where precision is lost  
void ppSumOther(xvector_t & v)
{
  xvector_t u( v.size(), 0);
  for (unsigned i = 0; i < v.size(); i++) {
    u[i] = 0;
    for (unsigned j = 0; j < v.size(); j++)
      if (i != j)
	u[i] += v[j];
  }
  v = u;
}


void transformByOptions(xvector_t & v, bool minusLogarithm, bool sumOther, string const & id)
{
  if (sumOther)
    ppSumOther(v);
  if (minusLogarithm)
    takeMinusLog(v, id);
}


// check identity of ids
void checkIds(string const & idVar, string const & idFac, unsigned lineCount)
{ 
  if (idVar != idFac)
    errorAbort("From main: problem with input number " + toString(lineCount) + ": idFac '" + idFac + "' and idVar '" + idVar + "' differ.");
}


// if facDataPtr not NULL, then read next line of potentials and reset dfg
void resetFactorPotential(FacData * facDataPtr, string const id, unsigned const lineCount, DFG & dfg)
{
  if (facDataPtr != NULL) {
    vector<xmatrix_t> facVec( facDataPtr->count() );
    string idFac;
    facDataPtr->next(idFac, facVec);
    checkIds(id, idFac, lineCount);
    dfg.resetFactorPotentials( facVec, facDataPtr->map() );
    dfg.consistencyCheck();
  }
}


// generate 2D vector of state symbols according to vector of stateMaps
vector< vector<string> > mkStateSymbolTable(vector<StateMapPtr_t> stateMapVec)
{ 
  vector< vector<string> > ssTable( stateMapVec.size() );
  for (unsigned i = 0; i < stateMapVec.size(); i++) {
    StateMap const & sm = *stateMapVec[i];
    for (unsigned j = 0; j < sm.stateCount(); j++)
      ssTable[i].push_back( sm.state2Symbol(j) );
  }
  return ssTable;
}
  

// convert vector of states to vector of symbols
void stateVecToSymbolVec(vector<StateMapPtr_t> const & stateMapVec, vector<state_t> const & maxVarStates, vector<symbol_t> & maxVarSymbols)
{
  assert( maxVarStates.size() == maxVarSymbols.size() );
  assert( maxVarStates.size() == stateMapVec.size() );
  for (unsigned i = 0; i < maxVarStates.size(); i++)
    maxVarSymbols[i] = stateMapVec[i]->state2Symbol( maxVarStates[i] );
}


// parse ppVarVecStr, which is of the form "X = a b c; Y = a b", and
// return variables (X and Y in this case) and states (a, b, and c in
// this case) as commonly indexed vectors.
void mkVarAndStateSymbolList(string const & varSpecStr,  vector<string> & varNames, vector< vector<symbol_t> > & varStates)
{
  varNames.clear();
  varStates.clear();
  vector<string> specs = split( strip(varSpecStr, "; \t"), ";");
  BOOST_FOREACH(string & s, specs) {
    // check for empty specs
    if (s.size() == 0 or strip(s).size() == 0)
      errorAbort("From mkVarAndStateSymbolList: Empty variable specification '" + s + "' found in this line:\n" + varSpecStr + "\n");

    vector<string> v = split(s, "=");
  
    //check
    if (not (v.size() == 1 or v.size() == 2) )
      errorAbort("From mkVarAndStateSymbolList: Error in specification of variable and state list.:\n" + varSpecStr + "\n");

    // add var name
    varNames.push_back( strip(v[0]) );

    // add var states (or none)
    vector<symbol_t> states;
    if (v.size() == 2) // states defined
      states = split( strip( v[1] ) );
    varStates.push_back(states);
  }
}


void writePostProbLegend(ostream & str, vector<string> const & varNames, vector< vector<symbol_t> > const & varStates)
{
  str << "#Definition of state order for each variable:" << endl;
  str << "#" << "NAME\t" << "ranVar" << "\tstate order ..." << endl;
  for (unsigned i = 0; i < varNames.size(); i++) {
    str << "#";
    writeNamedData(str, "name\t" + varNames[i], varStates[i]);
  }
}


int main(int argc, char * argv[])
{
  

  // option and argument variables
  string dfgSpecPrefix, stateMapsFile, factorPotentialsFile, variablesFile, factorGraphFile;  // specification files
  string varFile, facFile, postProbFile, normConstFile, maxProbStateFile; // input / ouput files
  string mpsVarVecStr, ppVarVecStr; // other options
  unsigned prec;
  bool minusLogarithm, ppSumOther;

  // SL: I added the following
  unsigned maxDepth;
  string ploidity;
  string model;

  // positional arguments (implemented as hidden options)
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("varFile", po::value<string>(& varFile), "Input variable file in named data format.")
    ("facFile", po::value<string>(& facFile), "Input factor file in named data format. Must use same identifiers in same order as varFile.");

  // define help message and options
  po::options_description visible(string("dfgEval allows implementation of discrete factor graphs and evaluates the probability of data sets under these models.\n\n")
				  + "  Usage: dfgEval [options] <inputVarData.tab> [inputFacData.tab]\n\n"
				  + "The arguments inputVarData.tab and inputFacData.tab are both in named data format.\n"
				  + "Allowed options");
  visible.add_options()
    ("help,h", "produce help message")
    ("ppFile,o", po::value<string>(& postProbFile)->default_value(""), "Calculate posterior probabilities for each state of each random variable and output to file.")
    ("ncFile,n", po::value<string>(& normConstFile)->default_value(""), "Calculate normalization constant output to file.")
    ("mpsFile,m", po::value<string>(& maxProbStateFile)->default_value(""), "Calculate most probable state for each random variable and output to file.")
    ("precision,p", po::value<unsigned>(& prec)->default_value(5), "Output precision of real numbers.")
    ("ppSumOther", po::bool_switch(& ppSumOther)->default_value(false), "For post probs, for each state output sum of post probs for all the other states for that variable. This retains precision for post probs very close to one.")
    ("minusLogarithm,l", po::bool_switch(& minusLogarithm)->default_value(false), "Output minus the natural logarithm of result values (program will terminate on negative results...).")
    ("mpsVars", po::value<string>(& mpsVarVecStr)->default_value(""), "Define the random variables for which the most probable state (mps) should be output. Default is to output the mps for all random variables. The specification string must be enclosed in citation marks and whitespace separated if it includes more than one random variable, e.g.: \"X Y\".")
    ("ppVars", po::value<string>(& ppVarVecStr)->default_value(""), "Define the random variables for which the posterior state probabilities (pp) should be calculated. Default is to output the pp for all states of all random variables (may generate much output!). Random variables are specified similar to mpsVars, but must be semicolon (';') separated. It is possible to only output pp's for certain states, in which case the following specification format is used: \"X=a b c; Y=a b\".")
    ("dfgSpecPrefix,s", po::value<string>(& dfgSpecPrefix)->default_value("./dfgSpec/"), "Prefix of DFG specification files..")
    ("factorGraphFile", po::value<string>(& factorGraphFile)->default_value("factorGraph.txt"), "Specification of the factor graph structure.")
    ("variablesFile", po::value<string>(& variablesFile)->default_value("variables.txt"), "Specification of the state map used by each variable.")
    ("stateMapFile", po::value<string>(& stateMapsFile)->default_value("stateMaps.txt"), "Specification of state maps.")
    ("facPotFile", po::value<string>(& factorPotentialsFile)->default_value("factorPotentials.txt"), "Specification of factor potentials.")
    // SL: I added the following
    ("maxDepth", po::value<unsigned>(& maxDepth)->default_value(200), "The maximum read depth. We expect all factorGraph.txt and variables.txt exist.")
    ("ploidity", po::value<string>(& ploidity)->default_value("diploid"), "The ploidity of the data.")
    ("model", po::value<string>(& model)->default_value("none"), "Specific model used (if any).");
  
  // SL: In the new version, we want to generate all DFGs for depth 1 to maxdepth
  // The files stateMapsFile and factorPotentialsFile depend on the ploidity parameter and the model used (if any).
  // diploid_stateMaps.txt and diploid_factorPotentials.txt or
  // haploid_stateMaps.txt and haploid_factorPotentials.txt
  // diploid_damage_stateMaps.txt and diploid_damage_factorPotentials.txt or
  // haploid_damage_stateMaps.txt and haploid_damage_factorPotentials.txt
  // The files factorGraphFile and variablesFile depend on the read depth and are called
  // depthN_factorGraph.txt and depthN_variables.txt
  // All files are placed in dfgSpecPrefix by default
  // We are going to generate the filenames on the fly and generate all the DFGs instead of reading only one.

  // setting up options parser
  po::options_description cmdline_options;
  cmdline_options.add(visible).add(hidden);

  po::positional_options_description p;
  p.add("varFile", 1).add("facFile", -1);

  po::variables_map vm;        
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);    

  // print help message
  if (vm.count("help")) {
    cout << visible << endl;
    return 1;
  }

  // check arguments
  if (vm.count("varFile") != 1)
    errorAbort("\nWrong number of arguments. Try -h for help");

  // set output precision at this point in case of xdouble type
#ifdef XNUMBER_IS_XDOUBLE
  xnumber_t::SetOutputPrecision(prec);
#endif

  //SL: I removed the following ...
  // read dfg
  /*  DfgInfo dfgInfo = readDfgInfo(dfgSpecPrefix + stateMapsFile, 
				dfgSpecPrefix + factorPotentialsFile, 
				dfgSpecPrefix + variablesFile, 
				dfgSpecPrefix + factorGraphFile);
  */

  // SL: ... and added this loop instead
  // Well, that didn't work. I'm doing it the hard and stupid way instead...
  // This should be dynamic ...  DfgInfo *dfgInfoList[maxDepth];
  DfgInfo *dfgInfoList[200];
  if(!model.empty()){
    model="_" + model;
  }
  string statemaps=dfgSpecPrefix + ploidity +  "_stateMaps.txt";
  string potentials=dfgSpecPrefix + ploidity + model + "_factorPotentials.txt";
  cerr<<statemaps<<endl<<potentials<<endl;
  FacData * facDataPtr = NULL;
  /*  if (facFile.size() != 0) {
    cout<<"facFile: "<<facFile<<endl;
    facDataPtr = new FacData(facFile, dfgInfo.facNames);
    }
  */

  // Create maxDepth DFGs and corresponding information. 
  // I do this the ugly way because loops and vectors didn't work. I should tidy this up at some point.
  // You can skip the next 400 lines...
  DfgInfo dfg1=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth1_variables.txt", dfgSpecPrefix + "depth1_factorGraph.txt");
  dfgInfoList[0]=&dfg1;
  DfgInfo dfg2=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth2_variables.txt", dfgSpecPrefix + "depth2_factorGraph.txt");
  dfgInfoList[1]=&dfg2;
  DfgInfo dfg3=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth3_variables.txt", dfgSpecPrefix + "depth3_factorGraph.txt");
  dfgInfoList[2]=&dfg3;
  DfgInfo dfg4=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth4_variables.txt", dfgSpecPrefix + "depth4_factorGraph.txt");
  dfgInfoList[3]=&dfg4;
  DfgInfo dfg5=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth5_variables.txt", dfgSpecPrefix + "depth5_factorGraph.txt");
  dfgInfoList[4]=&dfg5;
  DfgInfo dfg6=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth6_variables.txt", dfgSpecPrefix + "depth6_factorGraph.txt");
  dfgInfoList[5]=&dfg6;
  DfgInfo dfg7=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth7_variables.txt", dfgSpecPrefix + "depth7_factorGraph.txt");
  dfgInfoList[6]=&dfg7;
  DfgInfo dfg8=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth8_variables.txt", dfgSpecPrefix + "depth8_factorGraph.txt");
  dfgInfoList[7]=&dfg8;
  DfgInfo dfg9=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth9_variables.txt", dfgSpecPrefix + "depth9_factorGraph.txt");
  dfgInfoList[8]=&dfg9;
  DfgInfo dfg10=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth10_variables.txt", dfgSpecPrefix + "depth10_factorGraph.txt");
  dfgInfoList[9]=&dfg10;
  DfgInfo dfg11=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth11_variables.txt", dfgSpecPrefix + "depth11_factorGraph.txt");
  dfgInfoList[10]=&dfg11;
  DfgInfo dfg12=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth12_variables.txt", dfgSpecPrefix + "depth12_factorGraph.txt");
  dfgInfoList[11]=&dfg12;
  DfgInfo dfg13=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth13_variables.txt", dfgSpecPrefix + "depth13_factorGraph.txt");
  dfgInfoList[12]=&dfg13;
  DfgInfo dfg14=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth14_variables.txt", dfgSpecPrefix + "depth14_factorGraph.txt");
  dfgInfoList[13]=&dfg14;
  DfgInfo dfg15=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth15_variables.txt", dfgSpecPrefix + "depth15_factorGraph.txt");
  dfgInfoList[14]=&dfg15;
  DfgInfo dfg16=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth16_variables.txt", dfgSpecPrefix + "depth16_factorGraph.txt");
  dfgInfoList[15]=&dfg16;
  DfgInfo dfg17=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth17_variables.txt", dfgSpecPrefix + "depth17_factorGraph.txt");
  dfgInfoList[16]=&dfg17;
  DfgInfo dfg18=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth18_variables.txt", dfgSpecPrefix + "depth18_factorGraph.txt");
  dfgInfoList[17]=&dfg18;
  DfgInfo dfg19=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth19_variables.txt", dfgSpecPrefix + "depth19_factorGraph.txt");
  dfgInfoList[18]=&dfg19;
  DfgInfo dfg20=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth20_variables.txt", dfgSpecPrefix + "depth20_factorGraph.txt");
  dfgInfoList[19]=&dfg20;
  DfgInfo dfg21=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth21_variables.txt", dfgSpecPrefix + "depth21_factorGraph.txt");
  dfgInfoList[20]=&dfg21;
  DfgInfo dfg22=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth22_variables.txt", dfgSpecPrefix + "depth22_factorGraph.txt");
  dfgInfoList[21]=&dfg22;
  DfgInfo dfg23=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth23_variables.txt", dfgSpecPrefix + "depth23_factorGraph.txt");
  dfgInfoList[22]=&dfg23;
  DfgInfo dfg24=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth24_variables.txt", dfgSpecPrefix + "depth24_factorGraph.txt");
  dfgInfoList[23]=&dfg24;
  DfgInfo dfg25=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth25_variables.txt", dfgSpecPrefix + "depth25_factorGraph.txt");
  dfgInfoList[24]=&dfg25;
  DfgInfo dfg26=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth26_variables.txt", dfgSpecPrefix + "depth26_factorGraph.txt");
  dfgInfoList[25]=&dfg26;
  DfgInfo dfg27=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth27_variables.txt", dfgSpecPrefix + "depth27_factorGraph.txt");
  dfgInfoList[26]=&dfg27;
  DfgInfo dfg28=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth28_variables.txt", dfgSpecPrefix + "depth28_factorGraph.txt");
  dfgInfoList[27]=&dfg28;
  DfgInfo dfg29=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth29_variables.txt", dfgSpecPrefix + "depth29_factorGraph.txt");
  dfgInfoList[28]=&dfg29;
  DfgInfo dfg30=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth30_variables.txt", dfgSpecPrefix + "depth30_factorGraph.txt");
  dfgInfoList[29]=&dfg30;
  DfgInfo dfg31=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth31_variables.txt", dfgSpecPrefix + "depth31_factorGraph.txt");
  dfgInfoList[30]=&dfg31;
  DfgInfo dfg32=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth32_variables.txt", dfgSpecPrefix + "depth32_factorGraph.txt");
  dfgInfoList[31]=&dfg32;
  DfgInfo dfg33=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth33_variables.txt", dfgSpecPrefix + "depth33_factorGraph.txt");
  dfgInfoList[32]=&dfg33;
  DfgInfo dfg34=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth34_variables.txt", dfgSpecPrefix + "depth34_factorGraph.txt");
  dfgInfoList[33]=&dfg34;
  DfgInfo dfg35=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth35_variables.txt", dfgSpecPrefix + "depth35_factorGraph.txt");
  dfgInfoList[34]=&dfg35;
  DfgInfo dfg36=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth36_variables.txt", dfgSpecPrefix + "depth36_factorGraph.txt");
  dfgInfoList[35]=&dfg36;
  DfgInfo dfg37=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth37_variables.txt", dfgSpecPrefix + "depth37_factorGraph.txt");
  dfgInfoList[36]=&dfg37;
  DfgInfo dfg38=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth38_variables.txt", dfgSpecPrefix + "depth38_factorGraph.txt");
  dfgInfoList[37]=&dfg38;
  DfgInfo dfg39=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth39_variables.txt", dfgSpecPrefix + "depth39_factorGraph.txt");
  dfgInfoList[38]=&dfg39;
  DfgInfo dfg40=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth40_variables.txt", dfgSpecPrefix + "depth40_factorGraph.txt");
  dfgInfoList[39]=&dfg40;
  DfgInfo dfg41=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth41_variables.txt", dfgSpecPrefix + "depth41_factorGraph.txt");
  dfgInfoList[40]=&dfg41;
  DfgInfo dfg42=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth42_variables.txt", dfgSpecPrefix + "depth42_factorGraph.txt");
  dfgInfoList[41]=&dfg42;
  DfgInfo dfg43=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth43_variables.txt", dfgSpecPrefix + "depth43_factorGraph.txt");
  dfgInfoList[42]=&dfg43;
  DfgInfo dfg44=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth44_variables.txt", dfgSpecPrefix + "depth44_factorGraph.txt");
  dfgInfoList[43]=&dfg44;
  DfgInfo dfg45=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth45_variables.txt", dfgSpecPrefix + "depth45_factorGraph.txt");
  dfgInfoList[44]=&dfg45;
  DfgInfo dfg46=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth46_variables.txt", dfgSpecPrefix + "depth46_factorGraph.txt");
  dfgInfoList[45]=&dfg46;
  DfgInfo dfg47=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth47_variables.txt", dfgSpecPrefix + "depth47_factorGraph.txt");
  dfgInfoList[46]=&dfg47;
  DfgInfo dfg48=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth48_variables.txt", dfgSpecPrefix + "depth48_factorGraph.txt");
  dfgInfoList[47]=&dfg48;
  DfgInfo dfg49=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth49_variables.txt", dfgSpecPrefix + "depth49_factorGraph.txt");
  dfgInfoList[48]=&dfg49;
  DfgInfo dfg50=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth50_variables.txt", dfgSpecPrefix + "depth50_factorGraph.txt");
  dfgInfoList[49]=&dfg50;
  DfgInfo dfg51=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth51_variables.txt", dfgSpecPrefix + "depth51_factorGraph.txt");
  dfgInfoList[50]=&dfg51;
  DfgInfo dfg52=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth52_variables.txt", dfgSpecPrefix + "depth52_factorGraph.txt");
  dfgInfoList[51]=&dfg52;
  DfgInfo dfg53=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth53_variables.txt", dfgSpecPrefix + "depth53_factorGraph.txt");
  dfgInfoList[52]=&dfg53;
  DfgInfo dfg54=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth54_variables.txt", dfgSpecPrefix + "depth54_factorGraph.txt");
  dfgInfoList[53]=&dfg54;
  DfgInfo dfg55=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth55_variables.txt", dfgSpecPrefix + "depth55_factorGraph.txt");
  dfgInfoList[54]=&dfg55;
  DfgInfo dfg56=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth56_variables.txt", dfgSpecPrefix + "depth56_factorGraph.txt");
  dfgInfoList[55]=&dfg56;
  DfgInfo dfg57=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth57_variables.txt", dfgSpecPrefix + "depth57_factorGraph.txt");
  dfgInfoList[56]=&dfg57;
  DfgInfo dfg58=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth58_variables.txt", dfgSpecPrefix + "depth58_factorGraph.txt");
  dfgInfoList[57]=&dfg58;
  DfgInfo dfg59=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth59_variables.txt", dfgSpecPrefix + "depth59_factorGraph.txt");
  dfgInfoList[58]=&dfg59;
  DfgInfo dfg60=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth60_variables.txt", dfgSpecPrefix + "depth60_factorGraph.txt");
  dfgInfoList[59]=&dfg60;
  DfgInfo dfg61=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth61_variables.txt", dfgSpecPrefix + "depth61_factorGraph.txt");
  dfgInfoList[60]=&dfg61;
  DfgInfo dfg62=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth62_variables.txt", dfgSpecPrefix + "depth62_factorGraph.txt");
  dfgInfoList[61]=&dfg62;
  DfgInfo dfg63=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth63_variables.txt", dfgSpecPrefix + "depth63_factorGraph.txt");
  dfgInfoList[62]=&dfg63;
  DfgInfo dfg64=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth64_variables.txt", dfgSpecPrefix + "depth64_factorGraph.txt");
  dfgInfoList[63]=&dfg64;
  DfgInfo dfg65=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth65_variables.txt", dfgSpecPrefix + "depth65_factorGraph.txt");
  dfgInfoList[64]=&dfg65;
  DfgInfo dfg66=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth66_variables.txt", dfgSpecPrefix + "depth66_factorGraph.txt");
  dfgInfoList[65]=&dfg66;
  DfgInfo dfg67=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth67_variables.txt", dfgSpecPrefix + "depth67_factorGraph.txt");
  dfgInfoList[66]=&dfg67;
  DfgInfo dfg68=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth68_variables.txt", dfgSpecPrefix + "depth68_factorGraph.txt");
  dfgInfoList[67]=&dfg68;
  DfgInfo dfg69=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth69_variables.txt", dfgSpecPrefix + "depth69_factorGraph.txt");
  dfgInfoList[68]=&dfg69;
  DfgInfo dfg70=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth70_variables.txt", dfgSpecPrefix + "depth70_factorGraph.txt");
  dfgInfoList[69]=&dfg70;
  DfgInfo dfg71=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth71_variables.txt", dfgSpecPrefix + "depth71_factorGraph.txt");
  dfgInfoList[70]=&dfg71;
  DfgInfo dfg72=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth72_variables.txt", dfgSpecPrefix + "depth72_factorGraph.txt");
  dfgInfoList[71]=&dfg72;
  DfgInfo dfg73=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth73_variables.txt", dfgSpecPrefix + "depth73_factorGraph.txt");
  dfgInfoList[72]=&dfg73;
  DfgInfo dfg74=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth74_variables.txt", dfgSpecPrefix + "depth74_factorGraph.txt");
  dfgInfoList[73]=&dfg74;
  DfgInfo dfg75=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth75_variables.txt", dfgSpecPrefix + "depth75_factorGraph.txt");
  dfgInfoList[74]=&dfg75;
  DfgInfo dfg76=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth76_variables.txt", dfgSpecPrefix + "depth76_factorGraph.txt");
  dfgInfoList[75]=&dfg76;
  DfgInfo dfg77=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth77_variables.txt", dfgSpecPrefix + "depth77_factorGraph.txt");
  dfgInfoList[76]=&dfg77;
  DfgInfo dfg78=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth78_variables.txt", dfgSpecPrefix + "depth78_factorGraph.txt");
  dfgInfoList[77]=&dfg78;
  DfgInfo dfg79=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth79_variables.txt", dfgSpecPrefix + "depth79_factorGraph.txt");
  dfgInfoList[78]=&dfg79;
  DfgInfo dfg80=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth80_variables.txt", dfgSpecPrefix + "depth80_factorGraph.txt");
  dfgInfoList[79]=&dfg80;
  DfgInfo dfg81=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth81_variables.txt", dfgSpecPrefix + "depth81_factorGraph.txt");
  dfgInfoList[80]=&dfg81;
  DfgInfo dfg82=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth82_variables.txt", dfgSpecPrefix + "depth82_factorGraph.txt");
  dfgInfoList[81]=&dfg82;
  DfgInfo dfg83=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth83_variables.txt", dfgSpecPrefix + "depth83_factorGraph.txt");
  dfgInfoList[82]=&dfg83;
  DfgInfo dfg84=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth84_variables.txt", dfgSpecPrefix + "depth84_factorGraph.txt");
  dfgInfoList[83]=&dfg84;
  DfgInfo dfg85=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth85_variables.txt", dfgSpecPrefix + "depth85_factorGraph.txt");
  dfgInfoList[84]=&dfg85;
  DfgInfo dfg86=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth86_variables.txt", dfgSpecPrefix + "depth86_factorGraph.txt");
  dfgInfoList[85]=&dfg86;
  DfgInfo dfg87=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth87_variables.txt", dfgSpecPrefix + "depth87_factorGraph.txt");
  dfgInfoList[86]=&dfg87;
  DfgInfo dfg88=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth88_variables.txt", dfgSpecPrefix + "depth88_factorGraph.txt");
  dfgInfoList[87]=&dfg88;
  DfgInfo dfg89=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth89_variables.txt", dfgSpecPrefix + "depth89_factorGraph.txt");
  dfgInfoList[88]=&dfg89;
  DfgInfo dfg90=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth90_variables.txt", dfgSpecPrefix + "depth90_factorGraph.txt");
  dfgInfoList[89]=&dfg90;
  DfgInfo dfg91=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth91_variables.txt", dfgSpecPrefix + "depth91_factorGraph.txt");
  dfgInfoList[90]=&dfg91;
  DfgInfo dfg92=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth92_variables.txt", dfgSpecPrefix + "depth92_factorGraph.txt");
  dfgInfoList[91]=&dfg92;
  DfgInfo dfg93=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth93_variables.txt", dfgSpecPrefix + "depth93_factorGraph.txt");
  dfgInfoList[92]=&dfg93;
  DfgInfo dfg94=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth94_variables.txt", dfgSpecPrefix + "depth94_factorGraph.txt");
  dfgInfoList[93]=&dfg94;
  DfgInfo dfg95=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth95_variables.txt", dfgSpecPrefix + "depth95_factorGraph.txt");
  dfgInfoList[94]=&dfg95;
  DfgInfo dfg96=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth96_variables.txt", dfgSpecPrefix + "depth96_factorGraph.txt");
  dfgInfoList[95]=&dfg96;
  DfgInfo dfg97=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth97_variables.txt", dfgSpecPrefix + "depth97_factorGraph.txt");
  dfgInfoList[96]=&dfg97;
  DfgInfo dfg98=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth98_variables.txt", dfgSpecPrefix + "depth98_factorGraph.txt");
  dfgInfoList[97]=&dfg98;
  DfgInfo dfg99=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth99_variables.txt", dfgSpecPrefix + "depth99_factorGraph.txt");
  dfgInfoList[98]=&dfg99;
  DfgInfo dfg100=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth100_variables.txt", dfgSpecPrefix + "depth100_factorGraph.txt");
  dfgInfoList[99]=&dfg100;
  DfgInfo dfg101=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth101_variables.txt", dfgSpecPrefix + "depth101_factorGraph.txt");
  dfgInfoList[100]=&dfg101;
  DfgInfo dfg102=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth102_variables.txt", dfgSpecPrefix + "depth102_factorGraph.txt");
  dfgInfoList[101]=&dfg102;
  DfgInfo dfg103=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth103_variables.txt", dfgSpecPrefix + "depth103_factorGraph.txt");
  dfgInfoList[102]=&dfg103;
  DfgInfo dfg104=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth104_variables.txt", dfgSpecPrefix + "depth104_factorGraph.txt");
  dfgInfoList[103]=&dfg104;
  DfgInfo dfg105=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth105_variables.txt", dfgSpecPrefix + "depth105_factorGraph.txt");
  dfgInfoList[104]=&dfg105;
  DfgInfo dfg106=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth106_variables.txt", dfgSpecPrefix + "depth106_factorGraph.txt");
  dfgInfoList[105]=&dfg106;
  DfgInfo dfg107=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth107_variables.txt", dfgSpecPrefix + "depth107_factorGraph.txt");
  dfgInfoList[106]=&dfg107;
  DfgInfo dfg108=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth108_variables.txt", dfgSpecPrefix + "depth108_factorGraph.txt");
  dfgInfoList[107]=&dfg108;
  DfgInfo dfg109=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth109_variables.txt", dfgSpecPrefix + "depth109_factorGraph.txt");
  dfgInfoList[108]=&dfg109;
  DfgInfo dfg110=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth110_variables.txt", dfgSpecPrefix + "depth110_factorGraph.txt");
  dfgInfoList[109]=&dfg110;
  DfgInfo dfg111=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth111_variables.txt", dfgSpecPrefix + "depth111_factorGraph.txt");
  dfgInfoList[110]=&dfg111;
  DfgInfo dfg112=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth112_variables.txt", dfgSpecPrefix + "depth112_factorGraph.txt");
  dfgInfoList[111]=&dfg112;
  DfgInfo dfg113=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth113_variables.txt", dfgSpecPrefix + "depth113_factorGraph.txt");
  dfgInfoList[112]=&dfg113;
  DfgInfo dfg114=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth114_variables.txt", dfgSpecPrefix + "depth114_factorGraph.txt");
  dfgInfoList[113]=&dfg114;
  DfgInfo dfg115=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth115_variables.txt", dfgSpecPrefix + "depth115_factorGraph.txt");
  dfgInfoList[114]=&dfg115;
  DfgInfo dfg116=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth116_variables.txt", dfgSpecPrefix + "depth116_factorGraph.txt");
  dfgInfoList[115]=&dfg116;
  DfgInfo dfg117=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth117_variables.txt", dfgSpecPrefix + "depth117_factorGraph.txt");
  dfgInfoList[116]=&dfg117;
  DfgInfo dfg118=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth118_variables.txt", dfgSpecPrefix + "depth118_factorGraph.txt");
  dfgInfoList[117]=&dfg118;
  DfgInfo dfg119=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth119_variables.txt", dfgSpecPrefix + "depth119_factorGraph.txt");
  dfgInfoList[118]=&dfg119;
  DfgInfo dfg120=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth120_variables.txt", dfgSpecPrefix + "depth120_factorGraph.txt");
  dfgInfoList[119]=&dfg120;
  DfgInfo dfg121=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth121_variables.txt", dfgSpecPrefix + "depth121_factorGraph.txt");
  dfgInfoList[120]=&dfg121;
  DfgInfo dfg122=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth122_variables.txt", dfgSpecPrefix + "depth122_factorGraph.txt");
  dfgInfoList[121]=&dfg122;
  DfgInfo dfg123=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth123_variables.txt", dfgSpecPrefix + "depth123_factorGraph.txt");
  dfgInfoList[122]=&dfg123;
  DfgInfo dfg124=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth124_variables.txt", dfgSpecPrefix + "depth124_factorGraph.txt");
  dfgInfoList[123]=&dfg124;
  DfgInfo dfg125=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth125_variables.txt", dfgSpecPrefix + "depth125_factorGraph.txt");
  dfgInfoList[124]=&dfg125;
  DfgInfo dfg126=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth126_variables.txt", dfgSpecPrefix + "depth126_factorGraph.txt");
  dfgInfoList[125]=&dfg126;
  DfgInfo dfg127=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth127_variables.txt", dfgSpecPrefix + "depth127_factorGraph.txt");
  dfgInfoList[126]=&dfg127;
  DfgInfo dfg128=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth128_variables.txt", dfgSpecPrefix + "depth128_factorGraph.txt");
  dfgInfoList[127]=&dfg128;
  DfgInfo dfg129=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth129_variables.txt", dfgSpecPrefix + "depth129_factorGraph.txt");
  dfgInfoList[128]=&dfg129;
  DfgInfo dfg130=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth130_variables.txt", dfgSpecPrefix + "depth130_factorGraph.txt");
  dfgInfoList[129]=&dfg130;
  DfgInfo dfg131=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth131_variables.txt", dfgSpecPrefix + "depth131_factorGraph.txt");
  dfgInfoList[130]=&dfg131;
  DfgInfo dfg132=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth132_variables.txt", dfgSpecPrefix + "depth132_factorGraph.txt");
  dfgInfoList[131]=&dfg132;
  DfgInfo dfg133=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth133_variables.txt", dfgSpecPrefix + "depth133_factorGraph.txt");
  dfgInfoList[132]=&dfg133;
  DfgInfo dfg134=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth134_variables.txt", dfgSpecPrefix + "depth134_factorGraph.txt");
  dfgInfoList[133]=&dfg134;
  DfgInfo dfg135=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth135_variables.txt", dfgSpecPrefix + "depth135_factorGraph.txt");
  dfgInfoList[134]=&dfg135;
  DfgInfo dfg136=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth136_variables.txt", dfgSpecPrefix + "depth136_factorGraph.txt");
  dfgInfoList[135]=&dfg136;
  DfgInfo dfg137=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth137_variables.txt", dfgSpecPrefix + "depth137_factorGraph.txt");
  dfgInfoList[136]=&dfg137;
  DfgInfo dfg138=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth138_variables.txt", dfgSpecPrefix + "depth138_factorGraph.txt");
  dfgInfoList[137]=&dfg138;
  DfgInfo dfg139=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth139_variables.txt", dfgSpecPrefix + "depth139_factorGraph.txt");
  dfgInfoList[138]=&dfg139;
  DfgInfo dfg140=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth140_variables.txt", dfgSpecPrefix + "depth140_factorGraph.txt");
  dfgInfoList[139]=&dfg140;
  DfgInfo dfg141=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth141_variables.txt", dfgSpecPrefix + "depth141_factorGraph.txt");
  dfgInfoList[140]=&dfg141;
  DfgInfo dfg142=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth142_variables.txt", dfgSpecPrefix + "depth142_factorGraph.txt");
  dfgInfoList[141]=&dfg142;
  DfgInfo dfg143=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth143_variables.txt", dfgSpecPrefix + "depth143_factorGraph.txt");
  dfgInfoList[142]=&dfg143;
  DfgInfo dfg144=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth144_variables.txt", dfgSpecPrefix + "depth144_factorGraph.txt");
  dfgInfoList[143]=&dfg144;
  DfgInfo dfg145=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth145_variables.txt", dfgSpecPrefix + "depth145_factorGraph.txt");
  dfgInfoList[144]=&dfg145;
  DfgInfo dfg146=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth146_variables.txt", dfgSpecPrefix + "depth146_factorGraph.txt");
  dfgInfoList[145]=&dfg146;
  DfgInfo dfg147=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth147_variables.txt", dfgSpecPrefix + "depth147_factorGraph.txt");
  dfgInfoList[146]=&dfg147;
  DfgInfo dfg148=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth148_variables.txt", dfgSpecPrefix + "depth148_factorGraph.txt");
  dfgInfoList[147]=&dfg148;
  DfgInfo dfg149=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth149_variables.txt", dfgSpecPrefix + "depth149_factorGraph.txt");
  dfgInfoList[148]=&dfg149;
  DfgInfo dfg150=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth150_variables.txt", dfgSpecPrefix + "depth150_factorGraph.txt");
  dfgInfoList[149]=&dfg150;
  DfgInfo dfg151=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth151_variables.txt", dfgSpecPrefix + "depth151_factorGraph.txt");
  dfgInfoList[150]=&dfg151;
  DfgInfo dfg152=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth152_variables.txt", dfgSpecPrefix + "depth152_factorGraph.txt");
  dfgInfoList[151]=&dfg152;
  DfgInfo dfg153=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth153_variables.txt", dfgSpecPrefix + "depth153_factorGraph.txt");
  dfgInfoList[152]=&dfg153;
  DfgInfo dfg154=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth154_variables.txt", dfgSpecPrefix + "depth154_factorGraph.txt");
  dfgInfoList[153]=&dfg154;
  DfgInfo dfg155=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth155_variables.txt", dfgSpecPrefix + "depth155_factorGraph.txt");
  dfgInfoList[154]=&dfg155;
  DfgInfo dfg156=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth156_variables.txt", dfgSpecPrefix + "depth156_factorGraph.txt");
  dfgInfoList[155]=&dfg156;
  DfgInfo dfg157=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth157_variables.txt", dfgSpecPrefix + "depth157_factorGraph.txt");
  dfgInfoList[156]=&dfg157;
  DfgInfo dfg158=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth158_variables.txt", dfgSpecPrefix + "depth158_factorGraph.txt");
  dfgInfoList[157]=&dfg158;
  DfgInfo dfg159=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth159_variables.txt", dfgSpecPrefix + "depth159_factorGraph.txt");
  dfgInfoList[158]=&dfg159;
  DfgInfo dfg160=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth160_variables.txt", dfgSpecPrefix + "depth160_factorGraph.txt");
  dfgInfoList[159]=&dfg160;
  DfgInfo dfg161=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth161_variables.txt", dfgSpecPrefix + "depth161_factorGraph.txt");
  dfgInfoList[160]=&dfg161;
  DfgInfo dfg162=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth162_variables.txt", dfgSpecPrefix + "depth162_factorGraph.txt");
  dfgInfoList[161]=&dfg162;
  DfgInfo dfg163=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth163_variables.txt", dfgSpecPrefix + "depth163_factorGraph.txt");
  dfgInfoList[162]=&dfg163;
  DfgInfo dfg164=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth164_variables.txt", dfgSpecPrefix + "depth164_factorGraph.txt");
  dfgInfoList[163]=&dfg164;
  DfgInfo dfg165=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth165_variables.txt", dfgSpecPrefix + "depth165_factorGraph.txt");
  dfgInfoList[164]=&dfg165;
  DfgInfo dfg166=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth166_variables.txt", dfgSpecPrefix + "depth166_factorGraph.txt");
  dfgInfoList[165]=&dfg166;
  DfgInfo dfg167=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth167_variables.txt", dfgSpecPrefix + "depth167_factorGraph.txt");
  dfgInfoList[166]=&dfg167;
  DfgInfo dfg168=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth168_variables.txt", dfgSpecPrefix + "depth168_factorGraph.txt");
  dfgInfoList[167]=&dfg168;
  DfgInfo dfg169=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth169_variables.txt", dfgSpecPrefix + "depth169_factorGraph.txt");
  dfgInfoList[168]=&dfg169;
  DfgInfo dfg170=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth170_variables.txt", dfgSpecPrefix + "depth170_factorGraph.txt");
  dfgInfoList[169]=&dfg170;
  DfgInfo dfg171=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth171_variables.txt", dfgSpecPrefix + "depth171_factorGraph.txt");
  dfgInfoList[170]=&dfg171;
  DfgInfo dfg172=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth172_variables.txt", dfgSpecPrefix + "depth172_factorGraph.txt");
  dfgInfoList[171]=&dfg172;
  DfgInfo dfg173=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth173_variables.txt", dfgSpecPrefix + "depth173_factorGraph.txt");
  dfgInfoList[172]=&dfg173;
  DfgInfo dfg174=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth174_variables.txt", dfgSpecPrefix + "depth174_factorGraph.txt");
  dfgInfoList[173]=&dfg174;
  DfgInfo dfg175=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth175_variables.txt", dfgSpecPrefix + "depth175_factorGraph.txt");
  dfgInfoList[174]=&dfg175;
  DfgInfo dfg176=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth176_variables.txt", dfgSpecPrefix + "depth176_factorGraph.txt");
  dfgInfoList[175]=&dfg176;
  DfgInfo dfg177=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth177_variables.txt", dfgSpecPrefix + "depth177_factorGraph.txt");
  dfgInfoList[176]=&dfg177;
  DfgInfo dfg178=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth178_variables.txt", dfgSpecPrefix + "depth178_factorGraph.txt");
  dfgInfoList[177]=&dfg178;
  DfgInfo dfg179=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth179_variables.txt", dfgSpecPrefix + "depth179_factorGraph.txt");
  dfgInfoList[178]=&dfg179;
  DfgInfo dfg180=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth180_variables.txt", dfgSpecPrefix + "depth180_factorGraph.txt");
  dfgInfoList[179]=&dfg180;
  DfgInfo dfg181=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth181_variables.txt", dfgSpecPrefix + "depth181_factorGraph.txt");
  dfgInfoList[180]=&dfg181;
  DfgInfo dfg182=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth182_variables.txt", dfgSpecPrefix + "depth182_factorGraph.txt");
  dfgInfoList[181]=&dfg182;
  DfgInfo dfg183=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth183_variables.txt", dfgSpecPrefix + "depth183_factorGraph.txt");
  dfgInfoList[182]=&dfg183;
  DfgInfo dfg184=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth184_variables.txt", dfgSpecPrefix + "depth184_factorGraph.txt");
  dfgInfoList[183]=&dfg184;
  DfgInfo dfg185=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth185_variables.txt", dfgSpecPrefix + "depth185_factorGraph.txt");
  dfgInfoList[184]=&dfg185;
  DfgInfo dfg186=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth186_variables.txt", dfgSpecPrefix + "depth186_factorGraph.txt");
  dfgInfoList[185]=&dfg186;
  DfgInfo dfg187=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth187_variables.txt", dfgSpecPrefix + "depth187_factorGraph.txt");
  dfgInfoList[186]=&dfg187;
  DfgInfo dfg188=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth188_variables.txt", dfgSpecPrefix + "depth188_factorGraph.txt");
  dfgInfoList[187]=&dfg188;
  DfgInfo dfg189=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth189_variables.txt", dfgSpecPrefix + "depth189_factorGraph.txt");
  dfgInfoList[188]=&dfg189;
  DfgInfo dfg190=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth190_variables.txt", dfgSpecPrefix + "depth190_factorGraph.txt");
  dfgInfoList[189]=&dfg190;
  DfgInfo dfg191=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth191_variables.txt", dfgSpecPrefix + "depth191_factorGraph.txt");
  dfgInfoList[190]=&dfg191;
  DfgInfo dfg192=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth192_variables.txt", dfgSpecPrefix + "depth192_factorGraph.txt");
  dfgInfoList[191]=&dfg192;
  DfgInfo dfg193=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth193_variables.txt", dfgSpecPrefix + "depth193_factorGraph.txt");
  dfgInfoList[192]=&dfg193;
  DfgInfo dfg194=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth194_variables.txt", dfgSpecPrefix + "depth194_factorGraph.txt");
  dfgInfoList[193]=&dfg194;
  DfgInfo dfg195=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth195_variables.txt", dfgSpecPrefix + "depth195_factorGraph.txt");
  dfgInfoList[194]=&dfg195;
  DfgInfo dfg196=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth196_variables.txt", dfgSpecPrefix + "depth196_factorGraph.txt");
  dfgInfoList[195]=&dfg196;
  DfgInfo dfg197=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth197_variables.txt", dfgSpecPrefix + "depth197_factorGraph.txt");
  dfgInfoList[196]=&dfg197;
  DfgInfo dfg198=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth198_variables.txt", dfgSpecPrefix + "depth198_factorGraph.txt");
  dfgInfoList[197]=&dfg198;
  DfgInfo dfg199=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth199_variables.txt", dfgSpecPrefix + "depth199_factorGraph.txt");
  dfgInfoList[198]=&dfg199;
  DfgInfo dfg200=readDfgInfo(statemaps, potentials,dfgSpecPrefix + "depth200_variables.txt", dfgSpecPrefix + "depth200_factorGraph.txt");
  dfgInfoList[199]=&dfg200;
  // Done with DFGs...


  // open output stream - Fixed (SL)
  ofstream ppF;
  if (postProbFile != "-")
    openOutFile(ppF, postProbFile);
  ostream & ppStr = ( ppF.is_open() ) ? ppF : cout;

  // pp init data structures - Fixed (SL)
  vector<string>  ppVarNames;
  vector< vector<symbol_t> >  ppVarStates;
  mkVarAndStateSymbolList(ppVarVecStr,  ppVarNames, ppVarStates); 

  // We need a DFG to set some variables
  DfgInfo dfgInfo=*(dfgInfoList[maxDepth-1]);

  // pp output data structures
  vector<xvector_t> variableMarginals;
  initGenericVariableMarginals(variableMarginals, dfgInfo.dfg);

  // define which random variables and which states to output post probs for
  assert( ppVarNames.size() == ppVarStates.size() );
  vector<unsigned> ppVarMap = mkSubsetMap(dfgInfo.varNames, ppVarNames);

  // state maps
  vector< vector<string> > ssTable = mkStateSymbolTable(dfgInfo.stateMapVec); 
  vector< vector<unsigned> > ppVarStateMap;
  for (unsigned i = 0; i < ppVarStates.size(); i++) {
    if (ppVarStates[i].size() == 0)      ppVarStates[i] = ssTable[ ppVarMap[i] ];
    ppVarStateMap.push_back( mkSubsetMap( ssTable[ ppVarMap[i] ], ppVarStates[i] ) );
  }
  writeNamedData(ppStr, "NAME:\tranVar", ppVarStates[0]);

  // variables needed in data loop
  string idVar;
  vector<symbol_t> varVec;
  unsigned lineCount = 1;
  int depth;

  // Set the mapping using the largest DFG
  VarData varData(varFile, dfgInfo.varNames);
  vector<unsigned> theFullMap=varData.map();
  ifstream input(varFile.c_str());
  string myline;

  // Skip the first line with NAME: ...
  getline(input,myline);

  while ( getline(input, myline) ) {
    vector<string> elements=split(myline,'\t');
    idVar=elements[0];
    vector<symbol_t> varVec( elements.size()-1 );
    vector<unsigned> theMap( elements.size()-1 );
    
    // Find depth at the next position
    depth=min((int)maxDepth,atoi( (idVar.substr(idVar.find_last_of("_")+1)).c_str()) );
    for(int i=0;i<=depth;i++){
      varVec[i]=elements[i+1];
      theMap[i]=theFullMap[i];
    }

    DfgInfo dfgInfo=*(dfgInfoList[depth-1]);
    stateMaskVec_t stateMasks( dfgInfo.varNames.size() );

    dfgInfo.stateMaskMapSet.symbols2StateMasks(stateMasks, varVec, theMap);
    /*    for(int i=0;i<stateMasks.size();i++){
      cout<<stateMasks[i]<<" ";
    }
    cout<<endl;
    for(int i=0;i<varVec.size();i++){
      cout<<varVec[i]<<" ";
    }
    cout<<endl;

    for(int i=0;i<theMap.size();i++){
      cout<<theMap[i]<<" ";
    }
    cout<<endl;
    */

    dfgInfo.dfg.runSumProduct(stateMasks);  

    dfgInfo.dfg.calcVariableMarginals(variableMarginals, stateMasks);
    for (unsigned i = 0; i < ppVarNames.size(); i++) {
      xvector_t ppVec = variableMarginals[ ppVarMap[i] ];
      transformByOptions(ppVec, minusLogarithm, ppSumOther, idVar);
      writeNamedData(ppStr, idVar + "\t" + ppVarNames[i], mkSubset(toStdVector(ppVec), ppVarStateMap[i]), prec);
    }
    

    lineCount++;
  }

  // clean up
  if (facDataPtr != NULL)
    delete facDataPtr;
  input.close();

  return 0;
}

