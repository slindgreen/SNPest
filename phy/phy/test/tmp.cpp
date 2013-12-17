#include "phy/Grammar.h"

using namespace phy;

int main(void)
{
  Grammar g = readGrammar("data/grammarOneState.txt");

  // emit models
  xvector_t missingData(4, 1);
  vector<xvector_t> sngEmit(1, missingData);

  g.resetEmissions(sngEmit);

  cout << g.inside() << endl;
}
