#ifndef __NewickParser_h
#define __NewickParser_h

#include <string>
#include <vector>
#include "NewickNode.h"
#include "NewickParserBase.h"


class NewickParserBase; //forward declaration


/** 
 * NewickParser class can hold the information specified in the
 * newick format.
 */
class  NewickParser : public NewickParserBase{
public:
  NewickNode *makeNode(string nm,
		       double bl,
		       NewickNode *p,
		       vector<NewickNode *> *dl);
protected:
};

#endif  //__NewickParser_h
