#include<string.h>
#include<stdio.h>

extern void rhd_dq_module_MOD_rhd_dq_init( \
		char* /* parfile */, char* /* modelfile */, \
		int* /* nc_p */, int* /* imodel */, \
		int /* l_parfile */, int /* l_modelfile */);

extern void get_dq_(int* /* dq */);

#define rhd_dq_init_ __rhd_dq_module_MOD_rhd_dq_init

void rhd_dq_init(char* parfile, char* modelfile, \
		int nc_p, int imodel) {
	rhd_dq_init_(parfile, modelfile, &nc_p, &imodel, \
			strlen(parfile), strlen(modelfile));

	return;
}

void get_dq(int dq) {
	get_dq_(&dq);
	return;
}
