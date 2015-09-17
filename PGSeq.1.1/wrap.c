#include<Python.h>
#include <malloc.h>

double * calexpression(int IsoNO,int Rep);

//导出函数
PyObject *getNormconst(PyObject* self,PyObject* args)
{
	PyObject* tu;
	if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &tu))
   		 return NULL;
   	getNormc(tu);
   	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *get_genedata(PyObject* self,PyObject* args)
{
	char *GeneDataP;
	char *GeneMapP;
	char *GeneLenP;
	int IsoNo;
	int ExonNo;

	if(!PyArg_ParseTuple(args,"sssii",&GeneDataP,&GeneMapP,&GeneLenP,&IsoNo,&ExonNo))
	  return NULL;
	getgenedata(GeneDataP,GeneMapP,GeneLenP,IsoNo,ExonNo);
	Py_INCREF(Py_None);
	return Py_None;
}

PyObject *pgseq_calexpression(PyObject* self,PyObject* args) 
{
	int i,IsoNo,Rep;
	double *c;
	
	if(!PyArg_ParseTuple(args,"ii",&IsoNo,&Rep))
	  return NULL;
	  
	PyObject* plist=PyList_New(2*IsoNo*Rep+2*Rep);
	assert(PyList_Check(plist));
	
	c=calexpression(IsoNo,Rep);
	
	for(i=0;i<2*IsoNo*Rep+2*Rep;i++)
	{
		PyList_SetItem(plist,i,Py_BuildValue("d",c[i]));	
	}
	free(c);
	
	return plist; 
}

//方法列表
static PyMethodDef exampleMethods[]=
{
	{"getNormc",getNormconst,METH_VARARGS,"getNormconst"},
	{"getgenedata",get_genedata,METH_VARARGS,"get genedata"},
	{"calexpression",pgseq_calexpression,METH_VARARGS,"PGSeq calexpression"},
	{NULL,NULL}
};
//初始化函数
void initexample()
{
	PyObject* m;
	m=Py_InitModule("example",exampleMethods);
}
