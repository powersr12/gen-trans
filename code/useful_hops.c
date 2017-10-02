#include "useful_hops.h"

//Graphene TB - stores 1st, 2nd and 3rd NNTB hopping params for graphene
//last 2 params are placeholders for spin indices during transport calculations

int graphene_NNTB_max_neigh[] = 	{3, 9, 12};
double graphene_NNTB_lowdis[] = 	{0.56, 0.98, 1.13};
double graphene_NNTB_highdis[] = 	{0.59, 1.02, 1.17};
double graphene_NNTB_shifts[] = 	{0.0, 3*(-0.0740741), 0.0};
double graphene_NNTB_zmin[] = 		{0.0, 0.0, 0.0};
double graphene_NNTB_zmax[] = 		{0.01, 0.01, 0.01};
double _Complex graphene_NNTB_hops[] = 	{-1.0, -0.0740741, -0.0666667};

gen_hop_params graphene_NNTB_hop_params = {3,  graphene_NNTB_max_neigh, graphene_NNTB_hops, 0, 0, graphene_NNTB_lowdis,graphene_NNTB_highdis,graphene_NNTB_shifts,graphene_NNTB_zmin,  graphene_NNTB_zmax, 0, 0.0, NULL, NULL, 0, 0  };





//AB stacked Bilayer Graphene TB - stores inplane NN, out-of-plane direct interlayer, out-of-plane skew, 

//need to be careful! skew arent as easy as sep. dependence!

int BLG_NNTB_max_neigh[] = 	{3, 4, 7};
double BLG_NNTB_lowdis[] = 	{0.56, 0.0, 0.56};
double BLG_NNTB_highdis[] = 	{0.59, 0.1, 0.59};
double BLG_NNTB_shifts[] = 	{0.0, 0.0, 0.0};
double BLG_NNTB_zmin[] = 		{0.0, 13.3, 13.3};
double BLG_NNTB_zmax[] = 		{0.01, 13.5, 13.5};
double _Complex BLG_NNTB_hops[] = 	{-1.0, 0.12, 0.0};

gen_hop_params BLG_NNTB_hop_params = {3,  BLG_NNTB_max_neigh, BLG_NNTB_hops, 0, 0, BLG_NNTB_lowdis, BLG_NNTB_highdis, BLG_NNTB_shifts, BLG_NNTB_zmin,  BLG_NNTB_zmax, 0, 0.0, NULL, NULL, 0, 0  };




//Graphene TB with SOC - stores 1st NNTB hopping params for graphene, and 2NN intrinsic SOC hopping terms (& later others)
//last 2 params are placeholders for spin indices during transport calculations
//order of hoppings : NNTB, lambda_I_A, lambda_I_B, ...

int graphene_SOCTB_max_neigh[] = 	{3, 9, 9};
double graphene_SOCTB_lowdis[] = 	{0.56, 0.98, 0.98};
double graphene_SOCTB_highdis[] = 	{0.59, 1.02, 1.02};
double graphene_SOCTB_shifts[] = 	{0.0, 0.0, 0.0};
double graphene_SOCTB_zmin[] = 		{0.0, 0.0, 0.0};
double graphene_SOCTB_zmax[] = 		{0.01, 0.01, 0.01};
double _Complex graphene_SOCTB_hops[] = 	{-1.0, -0.01*I, -0.01*I};

gen_hop_params graphene_SOCTB_hop_params = {3,  graphene_SOCTB_max_neigh, graphene_SOCTB_hops, 0, 0, graphene_SOCTB_lowdis,graphene_SOCTB_highdis,graphene_SOCTB_shifts,graphene_SOCTB_zmin,  graphene_SOCTB_zmax, 0, 0.0, NULL, NULL, 0, 0  };