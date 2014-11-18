//*****************************************************************************
//  _____  __  __ _____         _____ ____  
// |  __ \|  \/  |  __ \       |_   _/ __ \
// | |__) | \  / | |  | |        | || |  | |
// |  ___/| |\/| | |  | |        | || |  | |
// | |    | |  | | |__| |       _| || |__| |
// |_|    |_|  |_|_____/ _____ |_____\____/ _wrapper
//                      |_____|
//
// Input/Output wrapper for the Porta Model Data (PMD)
//
//*****************************************************************************
//   C99
//   Flavio Calvo:               Geneva, Locarno
//   2014-10-10
//*****************************************************************************
//
//------*************----------------------------------------------------------
/*
 *
 * This is a C wrapper to the Fortran 90 I/O module for the Porta Model Data.
 * Only high-level routines are made available here:
 *
 * void pmd_init_header(struct pmd_header_type **pmd_header);
 * 	-> Allocates and initializes pmd_header with default values.
 *
 * void pmd_destroy_header(struct pmd_header_type **pmd_header);
 * 	-> Destroys pmd_header and frees memory.
 *
 * int pmd_openrd(int* unit, const char* file, \
 * 		struct pmd_header_type *pmd_header);
 * 	-> Opens file for reading. Handle to file is saved in unit, and
 * 	   PMD header is inserted inside pmd_header structure.
 *
 * int pmd_openwr(int* unit, const char* file, \
 * 		struct pmd_header_type *pmd_header);
 * 	-> Opens file for writing and writes pmd_header. Handle to file is
 * 	   saved in unit.
 *
 * int pmd_openap(int* unit, const char* file, \
 * 		struct pmd_header_type *pmd_header);
 * 	-> Opens file for appending new boxes. Header is read into pmd_header
 * 	   structure and handle is saved in unit.
 *
 * void pmd_init_box(struct pmd_box_type **p_box, \
 * 		struct pmd_header_type *pmd_header);
 * 	-> Allocates and initializes pmd_box according to pmd_header.
 *
 * void pmd_destroy_box(struct pmd_box_type **p_box);
 * 	-> Destroys pmd_box and frees memory
 *
 * int pmd_rd_box(const int unit, struct pmd_box_type *box);
 * 	-> Reads box from file designed by unit. Error code is returned if EOF.
 *
 * void pmd_append_box(const int unit, struct pmd_box_type *box);
 * 	-> Appends box to file designed by unit.
 *
 * void pmd_rd_module_hd(const int unit, int8_t* header, int sz);
 * 	-> Reads module header from file designed by unit.
 *
 * void pmd_wr_module_hd(const int unit, const int8_t* header, int sz);
 * 	-> Writes module header to file designed by unit.
 *
 * void pmd_close(const int unit);
 * 	-> Closes file designed by unit.
 *
 * How to link
 * -----------
 *  First step is to compule the Fortran module:
 *
 *  $ gfortran -c pmd_io_module.f90
 *
 *  This will produce an object file pmd_io_module.o. From this file one
 *  needs to call from here some functions belonging to the pypmd module.
 *
 *  The symbols to refer to these functions are prefixed (for the gfortran
 *  compiler) by __pypmd_MOD_. One can find this looking for the symbols
 *  inside pmd_io_module.o if an other compiler is used:
 *
 *  $ objdump -x pmd_io_module.o
 *
 *  Then one needs to modify the #define macros a few lines below.
 *
 *  Moreover, when Fortran passes strings of type character(len=*) to
 *  subroutines, there is an extra hidden integer argument passed after all
 *  visible arguments (at least with gfortran). It might happen with some
 *  compilers that the length is passed just after the string, and in such
 *  case this wrapper should be modified carefully to work.
 *
 *  One can then compiler the wrapper and the program using the wrapper:
 *
 *  $ gcc -c pmd_io_wrapper.c
 *  $ gcc -c my_program.c
 *
 *  And finally one can link with gcc:
 *
 *  $ gcc -o my_program pmd_io_module.o pmd_io_wrapper.o my_program.o -lgfortran
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define openrd __pypmd_MOD_openrd
#define openwr __pypmd_MOD_openwr
#define openap __pypmd_MOD_openap
#define rd_nbox __pypmd_MOD_rd_nbox
#define append_nbox __pypmd_MOD_append_nbox
#define rd_module_hd __pypmd_MOD_rd_module_hd
#define wr_module_hd __pypmd_MOD_wr_module_hd
#define close __pypmd_MOD_close

#define LEN_FILE 80
#define LEN_MAGIC_STR 8
#define LEN_XC 8192
#define LEN_ATOMIC_MODULE 1023
#define LEN_PMD_COMMENT 4096

extern void openrd(int *unit, int8_t *file, int8_t *magic_str, \
		int8_t *endianness, int8_t *int_sz, int8_t *db_sz, int \
		*pmd_version, int *localtime, int8_t *periodic, float \
		*domain_sz, float *domain_origin, int *dimensions, \
		float *xc1, float *xc2, float *xc3, int *n_radtheta, \
		int *n_radphi, int8_t *atomic_module, \
		int8_t *pmd_comment, \
		int *module_hd_sz, int *unused_var, \
		int *stat, int len_file, int len_magic_str, \
		int len_atomic_module, int len_pmd_comment);

extern void openwr(int *unit, int *ltime, int8_t *file, int8_t *magic_str, \
		int8_t *endianness, int8_t *int_sz, int8_t *db_sz, int \
		*pmd_version, int *localtime, int8_t *periodic, float \
		*domain_sz, float *domain_origin, int *dimensions, \
		float *xc1, float *xc2, float *xc3, int *n_radtheta, \
		int *n_radphi, int8_t *atomic_module, \
		int8_t *pmd_comment, \
		int *module_hd_sz, int *unused_var, \
		int *stat, int len_file, int len_magic_str, \
		int len_atomic_module, int len_pmd_comment);

extern void openap(int *unit, int8_t *file, int8_t *magic_str, \
		int8_t *endianness, int8_t *int_sz, int8_t *db_sz, int \
		*pmd_version, int *localtime, int8_t *periodic, float \
		*domain_sz, float *domain_origin, int *dimensions, \
		float *xc1, float *xc2, float *xc3, int *n_radtheta, \
		int *n_radphi, int8_t *atomic_module, \
		int8_t *pmd_comment, \
		int *module_hd_sz, int *unused_var, \
		int *stat, int len_file, int len_magic_str, \
		int len_atomic_module, int len_pmd_comment);

extern void rd_nbox(const int *unit, double *box, int *iostat, \
		int *nx, int *ny, int *nz);

extern void append_nbox(const int *unit, double *box, \
		int *nx, int *ny, int *nz);

extern void rd_module_hd(const int *unit, int8_t *module_header, int *sz);

extern void wr_module_hd(const int *unit, const int8_t *module_header, \
		int len_module_header);

extern void close(const int *unit);

struct pmd_header_type {
	int8_t *magic_str;
	int8_t *endianness;
	int8_t *int_sz, *db_sz;
	int *pmd_version;
	int *localtime;
	int8_t *periodic;
	float *domain_sz, *domain_origin;
	int *dimensions;
	float *xc1, *xc2, *xc3;
	int *n_radtheta, *n_radphi;
	int8_t *atomic_module;
	int8_t *pmd_comment;
	int *module_hd_sz;
	int *unused_var;
};

struct pmd_box_type {
	int nx;
	int ny;
	int nz;
	double *box;
};

void pmd_init_header(struct pmd_header_type **pmd_header);
void pmd_destroy_header(struct pmd_header_type **pmd_header);
void pmd_setstr(int8_t *str, const char *val, int len);
int pmd_openrd(int* unit, const char* file, \
		struct pmd_header_type *pmd_header);
int pmd_openwr(int* unit, const char* file, \
		struct pmd_header_type *pmd_header);
int pmd_openap(int* unit, const char* file, \
		struct pmd_header_type *pmd_header);
void pmd_init_box(struct pmd_box_type **p_box, \
		struct pmd_header_type *pmd_header);
void pmd_destroy_box(struct pmd_box_type **p_box);
int pmd_rd_box(const int unit, struct pmd_box_type *box);
void pmd_append_box(const int unit, struct pmd_box_type *box);
void pmd_rd_module_hd(const int unit, int8_t* header, int sz);
void pmd_wr_module_hd(const int unit, const int8_t* header, int sz);
void pmd_close(const int unit);

void pmd_init_header(struct pmd_header_type **p_pmd_header){
	struct pmd_header_type *pmd_header;
	*p_pmd_header=malloc(sizeof(struct pmd_header_type));
	pmd_header=*p_pmd_header;
	pmd_header->magic_str=malloc(LEN_MAGIC_STR*sizeof(int8_t));
	pmd_setstr(pmd_header->magic_str, "portapmd", 0);
	pmd_header->endianness=malloc(sizeof(int8_t));
	*(pmd_header->endianness)=0;
	pmd_header->int_sz=malloc(sizeof(int8_t));
	*(pmd_header->int_sz)=4;
	pmd_header->db_sz=malloc(sizeof(int8_t));
	*(pmd_header->db_sz)=8;
	pmd_header->pmd_version=malloc(sizeof(int));
	*(pmd_header->pmd_version)=2;
	pmd_header->localtime=malloc(6*sizeof(int));
	memset(pmd_header->localtime, 0, 6*sizeof(int));
	pmd_header->periodic=malloc(2*sizeof(int8_t));
	memset(pmd_header->periodic, 0, 2*sizeof(int8_t));
	pmd_header->domain_sz=malloc(3*sizeof(float));
	memset(pmd_header->domain_sz, 0, 3*sizeof(float));
	pmd_header->domain_origin=malloc(3*sizeof(float));
	memset(pmd_header->domain_origin, 0, 3*sizeof(float));
	pmd_header->dimensions=malloc(3*sizeof(int));
	memset(pmd_header->dimensions, 0, 3*sizeof(int));
	pmd_header->xc1=malloc(LEN_XC*sizeof(float));
	memset(pmd_header->xc1, 0, LEN_XC*sizeof(float));
	pmd_header->xc2=malloc(LEN_XC*sizeof(float));
	memset(pmd_header->xc2, 0, LEN_XC*sizeof(float));
	pmd_header->xc3=malloc(LEN_XC*sizeof(float));
	memset(pmd_header->xc3, 0, LEN_XC*sizeof(float));
	pmd_header->n_radtheta=malloc(sizeof(int));
	*(pmd_header->n_radtheta)=0;
	pmd_header->n_radphi=malloc(sizeof(int));
	*(pmd_header->n_radphi)=0;
	pmd_header->atomic_module=malloc(LEN_ATOMIC_MODULE*sizeof(int8_t));
	memset(pmd_header->atomic_module, (const int)' ', LEN_ATOMIC_MODULE*sizeof(int8_t));
	pmd_header->pmd_comment=malloc(LEN_PMD_COMMENT*sizeof(int8_t));
	memset(pmd_header->pmd_comment, (const int)' ', LEN_PMD_COMMENT*sizeof(int8_t));
	pmd_header->module_hd_sz=malloc(sizeof(int));
	*(pmd_header->module_hd_sz)=0;
	pmd_header->unused_var=malloc(sizeof(int));
	*(pmd_header->unused_var)=0;
}

void pmd_destroy_header(struct pmd_header_type **p_pmd_header){
	struct pmd_header_type *pmd_header;
	pmd_header=*p_pmd_header;
	free(pmd_header->magic_str);
	free(pmd_header->endianness);
	free(pmd_header->int_sz);
	free(pmd_header->db_sz);
	free(pmd_header->pmd_version);
	free(pmd_header->localtime);
	free(pmd_header->periodic);
	free(pmd_header->domain_sz);
	free(pmd_header->domain_origin);
	free(pmd_header->dimensions);
	free(pmd_header->xc1);
	free(pmd_header->xc2);
	free(pmd_header->xc3);
	free(pmd_header->n_radtheta);
	free(pmd_header->n_radphi);
	free(pmd_header->atomic_module);
	free(pmd_header->pmd_comment);
	free(pmd_header->module_hd_sz);
	free(pmd_header->unused_var);
	free(pmd_header);
	*p_pmd_header=NULL;
}

void pmd_setstr(int8_t *str, const char *val, int len){
	int i;
	int end=0;

	for(i=0;;i++){
		if(val[i] == '\0') end=1;
		if(!end && (i<len || len==0)) str[i]=val[i];
		else if(i<len) str[i]=' ';
		else break;
	}
	return;
}

int pmd_openrd(int* unit, const char* file, \
		struct pmd_header_type *pmd_header){
	int stat;
	openrd(unit, (int8_t*) file, pmd_header->magic_str, \
			pmd_header->endianness, pmd_header->int_sz, \
			pmd_header->db_sz, pmd_header->pmd_version, \
			pmd_header->localtime, pmd_header->periodic, \
			pmd_header->domain_sz, pmd_header->domain_origin, \
			pmd_header->dimensions, pmd_header->xc1, \
			pmd_header->xc2, pmd_header->xc3, \
			pmd_header->n_radtheta, pmd_header->n_radphi, \
			pmd_header->atomic_module, pmd_header->pmd_comment, \
			pmd_header->module_hd_sz, pmd_header->unused_var, \
			&stat, LEN_FILE, LEN_MAGIC_STR, \
			LEN_ATOMIC_MODULE, LEN_PMD_COMMENT);
	return stat;
}

int pmd_openwr(int* unit, const char* file, \
		struct pmd_header_type *pmd_header){
	int *ltime, stat;
	ltime=malloc(6*sizeof(int));
	openwr(unit, ltime, (int8_t*) file, pmd_header->magic_str, \
			pmd_header->endianness, pmd_header->int_sz, \
			pmd_header->db_sz, pmd_header->pmd_version, \
			pmd_header->localtime, pmd_header->periodic, \
			pmd_header->domain_sz, pmd_header->domain_origin, \
			pmd_header->dimensions, pmd_header->xc1, \
			pmd_header->xc2, pmd_header->xc3, \
			pmd_header->n_radtheta, pmd_header->n_radphi, \
			pmd_header->atomic_module, pmd_header->pmd_comment, \
			pmd_header->module_hd_sz, pmd_header->unused_var, \
			&stat, LEN_FILE, LEN_MAGIC_STR, \
			LEN_ATOMIC_MODULE, LEN_PMD_COMMENT);
	return stat;
}

int pmd_openap(int* unit, const char* file, \
		struct pmd_header_type *pmd_header){
	int stat;
	openap(unit, (int8_t*) file, pmd_header->magic_str, \
			pmd_header->endianness, pmd_header->int_sz, \
			pmd_header->db_sz, pmd_header->pmd_version, \
			pmd_header->localtime, pmd_header->periodic, \
			pmd_header->domain_sz, pmd_header->domain_origin, \
			pmd_header->dimensions, pmd_header->xc1, \
			pmd_header->xc2, pmd_header->xc3, \
			pmd_header->n_radtheta, pmd_header->n_radphi, \
			pmd_header->atomic_module, pmd_header->pmd_comment, \
			pmd_header->module_hd_sz, pmd_header->unused_var, \
			&stat, LEN_FILE, LEN_MAGIC_STR, \
			LEN_ATOMIC_MODULE, LEN_PMD_COMMENT);
	return stat;
}

void pmd_init_box(struct pmd_box_type **p_box, \
		struct pmd_header_type *pmd_header){
	struct pmd_box_type *box;
	*p_box=malloc(sizeof(struct pmd_box_type));
	box=*p_box;
	box->nx=pmd_header->dimensions[0];
	box->ny=pmd_header->dimensions[1];
	box->nz=pmd_header->dimensions[2];
	box->box=malloc(box->nx*box->ny*box->nz*sizeof(double));
}

void pmd_destroy_box(struct pmd_box_type **p_box){
	struct pmd_box_type *box;
	box=*p_box;
	free(box->box);
	free(box);
	*p_box=NULL;
}

int pmd_rd_box(const int unit, struct pmd_box_type *box){
	int iostat;
	rd_nbox(&unit, box->box, &iostat, &box->nx, &box->ny, &box->nz);
	return iostat;
}

void pmd_append_box(const int unit, struct pmd_box_type *box){
	append_nbox(&unit, box->box, &box->nx, &box->ny, &box->nz);
}

void pmd_rd_module_hd(const int unit, int8_t* header, int sz){
	rd_module_hd(&unit, header, &sz);
}

void pmd_wr_module_hd(const int unit, const int8_t* header, int sz){
	wr_module_hd(&unit, header, sz);
}

void pmd_close(const int unit){
	close(&unit);
}

int main(int argc, char *argv[]){
	struct pmd_header_type *pmd_header;
	struct pmd_box_type *box, *box2;
	int unit;
	int8_t *module_header;
	int i,j,k;
	double arr[24];
	char file[]="testd.pmd";

	pmd_init_header(&pmd_header);

	// Make header
	printf("Make header...\n");
	*(pmd_header->endianness)=0;
	*(pmd_header->int_sz)=4;
	*(pmd_header->db_sz)=8;
	pmd_header->periodic[0]=1;
	pmd_header->periodic[1]=0;
	pmd_header->domain_sz[0]=20.;
	pmd_header->domain_sz[1]=30.;
	pmd_header->domain_sz[2]=40.;
	pmd_header->domain_origin[0]=1.;
	pmd_header->domain_origin[1]=2.;
	pmd_header->domain_origin[2]=3.;
	pmd_header->dimensions[0]=2.;
	pmd_header->dimensions[1]=3.;
	pmd_header->dimensions[2]=4.;
	pmd_header->xc1[0]=6.;
	pmd_header->xc1[1]=16.;
	pmd_header->xc2[0]=7.;
	pmd_header->xc2[1]=17.;
	pmd_header->xc2[2]=27.;
	pmd_header->xc3[0]=8.;
	pmd_header->xc3[1]=18.;
	pmd_header->xc3[2]=28.;
	pmd_header->xc3[3]=38.;
	*(pmd_header->n_radtheta)=2;
	*(pmd_header->n_radphi)=4;
	pmd_setstr(pmd_header->atomic_module, "PORTA Continuum", 0);
	pmd_setstr(pmd_header->pmd_comment, "The typical small box for testing", 0);
	module_header=malloc(847*sizeof(int8_t));
	pmd_setstr(module_header, "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed non risus. Suspendisse lectus tortor, dignissim sit amet, adipiscing nec, ultricies sed, dolor. Cras elementum ultrices diam. Maecenas ligula massa, varius a, semper congue, euismod non, mi. Proin porttitor, orci nec nonummy molestie, enim est eleifend mi, non fermentum diam nisl sit amet erat. Duis semper. Duis arcu massa, scelerisque vitae, consequat in, pretium a, enim. Pellentesque congue. Ut in risus volutpat libero pharetra tempor. Cras vestibulum bibendum augue. Praesent egestas leo in pede. Praesent blandit odio eu enim. Pellentesque sed dui ut augue blandit sodales. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia Curae; Aliquam nibh. Mauris ac mauris sed pede pellentesque fermentum. Maecenas adipiscing ante non diam sodales hendrerit.", 0);
	*(pmd_header->module_hd_sz)=847;

	// Make the first box
	printf("Make the first box...\n");
	pmd_init_box(&box, pmd_header);
	for(i=0;i<2;i++)
	for(j=0;j<3;j++)
	for(k=0;k<4;k++) box->box[i+2*j+6*k]=(double)(i+2*j+6*k);

	// Write PMD file
	printf("Write PMD file...\n");
	pmd_openwr(&unit, file, pmd_header);
	printf("Unit: %d.\n", unit);
	pmd_wr_module_hd(unit, module_header, *pmd_header->module_hd_sz);
	pmd_append_box(unit, box);
	// Inverse the box
	for(i=0;i<24;i++) box->box[i]=-box->box[i];
	pmd_append_box(unit, box);
	pmd_close(unit);

	// Open file and read
	pmd_destroy_header(&pmd_header);
	pmd_init_header(&pmd_header);
	pmd_openrd(&unit, file, pmd_header);
	module_header[0]=0;
	pmd_rd_module_hd(unit, module_header, *pmd_header->module_hd_sz);
	printf("%s\n", module_header);
	pmd_rd_box(unit, box);
	pmd_init_box(&box2, pmd_header);
	pmd_rd_box(unit, box2);
	pmd_close(unit);

	// Write file again
	strcpy(file, "teste.pmd");
	pmd_openwr(&unit, file, pmd_header);
	pmd_wr_module_hd(unit, module_header, *pmd_header->module_hd_sz);
	pmd_append_box(unit, box);
	pmd_append_box(unit, box2);
	pmd_close(unit);

	printf("PMD version: %d.\n", *(pmd_header->pmd_version));
	pmd_destroy_header(&pmd_header);

	return 0;
}
