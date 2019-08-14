#include "devices.h"
#include "connect.h"
#include "transport.h"
#include "useful_hops.h"
#include <stddef.h>


//Graphene TB - stores 1st, 2nd and 3rd NNTB hopping params for graphene

int graphene_NNTB_max_neigh[] = 	{3, 9, 12};
double graphene_NNTB_lowdis[] = 	{0.56, 0.98, 1.13};
double graphene_NNTB_highdis[] = 	{0.63, 1.02, 1.17};
double graphene_NNTB_shifts[] = 	{0.0, 3*(-0.0740741), 0.0};
double graphene_NNTB_zmin[] = 		{0.0, 0.0, 0.0};
double graphene_NNTB_zmax[] = 		{0.01, 0.01, 0.01};
double _Complex graphene_NNTB_hops[] = 	{-1.0, -0.0740741, -0.0666667};

gen_hop_params graphene_NNTB_hop_params = {3,  graphene_NNTB_max_neigh, graphene_NNTB_hops, 0, 0, graphene_NNTB_lowdis,graphene_NNTB_highdis,graphene_NNTB_shifts,graphene_NNTB_zmin,  graphene_NNTB_zmax, 0, 0.0, NULL, NULL  };





//AB stacked Bilayer Graphene TB - stores inplane NN, out-of-plane direct interlayer, out-of-plane skew, 

//need to be careful! skew arent as easy as sep. dependence!

int BLG_NNTB_max_neigh[] = 	{3, 4, 7};
double BLG_NNTB_lowdis[] = 	{0.56, 0.0, 0.56};
double BLG_NNTB_highdis[] = 	{0.59, 0.1, 0.59};
double BLG_NNTB_shifts[] = 	{0.0, 0.0, 0.0};
double BLG_NNTB_zmin[] = 		{0.0, 1.36, 1.36};
double BLG_NNTB_zmax[] = 		{0.01, 1.37, 1.37};
double _Complex BLG_NNTB_hops[] = 	{-1.0, 0.12, 0.0};

gen_hop_params BLG_NNTB_hop_params = {3,  BLG_NNTB_max_neigh, BLG_NNTB_hops, 0, 0, BLG_NNTB_lowdis, BLG_NNTB_highdis, BLG_NNTB_shifts, BLG_NNTB_zmin,  BLG_NNTB_zmax, 0, 0.0, NULL, NULL  };


//Multilayer system  
//Based on a certain number of connections:
//nearest neighbour ONLY within layer
//vertical, and first non-vertical interlayer (or equiv.) between layers
//associated hopping different to others (uses same formula for in and out of plane)
//these params only really used for connection and thresholds.... (probably)

int MLG_NNTB_max_neigh[] = 	{3, 9};
double MLG_NNTB_lowdis[] = 	{0.56, 0.0};
double MLG_NNTB_highdis[] = 	{0.59, 0.59};
double MLG_NNTB_shifts[] = 	{0.0, 0.0};
double MLG_NNTB_zmin[] = 		{0.0, 1.36};
double MLG_NNTB_zmax[] = 		{0.01, 1.37};
double _Complex MLG_NNTB_hops[] = 	{-1.0, 0.178};

gen_hop_params MLG_NNTB_hop_params = {2,  MLG_NNTB_max_neigh, MLG_NNTB_hops, 0, 0, MLG_NNTB_lowdis, MLG_NNTB_highdis, MLG_NNTB_shifts, MLG_NNTB_zmin,  MLG_NNTB_zmax, 0, 0.0, NULL, NULL  };