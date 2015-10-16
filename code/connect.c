#include "connect.h"

void device_connectivity (RectRedux *DeviceCell, cnxRulesFn *rule, void *rule_params, cnxProfile *cnxp)
{
    int N = *(DeviceCell->Ntot);
    int i, j, counter;
    
    (cnxp->site_cnxnum) = createIntArray(N);
    
    (cnxp->site_cnx) = createNonSquareIntMatrix(N, cnxp->max_neigh);
    
    for(i=0; i<N; i++)
    {
      counter=0;
      if((DeviceCell->siteinfo)[i][0] == 0)
      {
	for(j=0; j<N; j++)
	{
	  if((DeviceCell->siteinfo)[j][0] == 0)
	  {
	    if((*rule)(DeviceCell, rule_params, i, j) == 0)
	    {
	      (cnxp->site_cnx)[i][counter] = j;
	      counter++;
	    }
	  }
	}
      }
      (cnxp->site_cnxnum)[i] = counter;
    }
    
}

void printConnectivity (RectRedux *DeviceCell, cnxProfile *cnxp)
{
  int i, j;
  int N = *(DeviceCell->Ntot);

  for(i=0; i<N; i++)
  {
    printf("%lf	%lf\n", (DeviceCell->pos)[i][0], (DeviceCell->pos)[i][1]);
    
    for(j=0; j<(cnxp->site_cnxnum)[i]; j++)
    {
      printf("%lf	%lf\n%lf	%lf\n",(DeviceCell->pos)[cnxp->site_cnx[i][j]][0], (DeviceCell->pos)[cnxp->site_cnx[i][j]][1], (DeviceCell->pos)[i][0], (DeviceCell->pos)[i][1]  );
    }
    
    printf("\n");
  }
  
}


//generates a starting cell for the cell splitter algorithm
//config=0 is the standard 2 nanoribbon lead case.
//starting cell means the cell whose GFs we want at the end of the recursion
//depending on number of sweeps it can be the first or last actual cell
//for a simple transport sweep, it is the cell connecting to the final (RHS) lead
//also generates some of the arrays in cellDivision structure
void genStartingCell (RectRedux *DeviceCell, cellDivision *cellinfo, int config, void *start_params)
{
    //cellinfo->num_cells = 2;
    int i, j, k, l, m, n;
    int dim = *(DeviceCell->Nrem);
    int alldim =2*(DeviceCell->length)*(DeviceCell->length2);
    double xmin, xmax;
    //printf("%d	%d\n", dim, alldim);
    //cellDivision.cells_site_order =
    
    cellinfo->cells_site_order = createIntArray(dim);
    cellinfo->sites_by_cell = createIntArray(alldim);

    //initialise cells_site_order elements to -1
    for(i=0; i<dim; i++)
    {
      (cellinfo->cells_site_order)[i] = -1;
    }
    
    //initialise sites_by_cell
    //-2 = not in lattice, -1 in lattice, not yet assigned
    for(i=0; i<alldim; i++)
    {
      if( (DeviceCell->siteinfo)[i][0] == 1)
	(cellinfo->sites_by_cell)[i] = -2;
      
      if( (DeviceCell->siteinfo)[i][0] == 0)
	(cellinfo->sites_by_cell)[i] = -1;
    }
    
    
    
    if(config==0)	// standard L/R leads same width as device. 
			// cell 0 is the cell connecting to the RHS probe
			// the "group" is the LHS cell
    {
      j=0;
      for(i=0; i< 2*(DeviceCell->length); i++)
      {
	if( (DeviceCell->siteinfo)[alldim - 2*(DeviceCell->length) + i][0] == 0)
	{
	    (cellinfo->cells_site_order)[j] = alldim - 2*(DeviceCell->length) + i;
	    j++;
	    
	    (cellinfo->sites_by_cell)[alldim - 2*(DeviceCell->length) + i] = 0;
	}
	
		
      }
      
      cellinfo->cell1dim = j;
      
      //place the sites that will connect to the other lead in a group together
      cellinfo->group_dim = 2*(DeviceCell->length);
      cellinfo->group_sites = createIntArray(cellinfo->group_dim);
      
      for(i=0; i< 2*(DeviceCell->length); i++)
      {
	cellinfo->group_sites[i]=i;
      }
      
    }
    
    else if(config==1) //just the sites that actually connect to leads, and not the full conventional cells
    {
      j=0;
      
      if( (DeviceCell->geo) == 0)	//zigzag case
      {
	for(i=0; i< 2*(DeviceCell->length); i++)
	{
	  m= (i % 4); //which atom in 4 atom unit cell is i
	  
	  if(m==0 || m==3)
	  {
	 
	    if( (DeviceCell->siteinfo)[alldim - 2*(DeviceCell->length) + i][0] == 0)
	    {
		(cellinfo->cells_site_order)[j] = alldim - 2*(DeviceCell->length) + i;
		j++;
		
		(cellinfo->sites_by_cell)[alldim - 2*(DeviceCell->length) + i] = 0;
	    }
	
	  }
	}
	cellinfo->cell1dim = j;  
	//printf("#cell1 %d\n", cellinfo->cell1dim);
	
	//group sites
	l=0;
	for(i=0; i< 2*(DeviceCell->length); i++)
	{
	  m= (i % 4); //which atom in 4 atom unit cell is i
	  
	  if(m==1 || m==2)
	  {
	    if( (DeviceCell->siteinfo)[i][0] == 0)
	    {
	     l++; 
	    }
	  }
	}
	cellinfo->group_dim = l;
	//printf("#group %d\n", cellinfo->group_dim);

	cellinfo->group_sites = createIntArray(cellinfo->group_dim);
	l=0;
	for(i=0; i< 2*(DeviceCell->length); i++)
	{
	  m= (i % 4); //which atom in 4 atom unit cell is i
	  
	  if(m==1 || m==2)
	  {
	    if( (DeviceCell->siteinfo)[i][0] == 0)
	    {
	      (cellinfo->group_sites)[l] = i;
	      l++;
	    }
	  }
	}
	
      
      }
      
      
      if( (DeviceCell->geo) == 1)	//armchair case
      {
	for(i=0; i< (DeviceCell->length); i++)
	{
	  m= (i % 2); //which atom in 2 atom unit cell is i
	  
	  if(m==0 )
	  {
	 
	    if( (DeviceCell->siteinfo)[alldim - (DeviceCell->length) + i][0] == 0)
	    {
		(cellinfo->cells_site_order)[j] = alldim -(DeviceCell->length) + i;
		j++;
		
		(cellinfo->sites_by_cell)[alldim - (DeviceCell->length) + i] = 0;
	    }
	
	  }
	}
	cellinfo->cell1dim = j;  
	//printf("#cell1 %d\n", cellinfo->cell1dim);
	
	//group sites
	l=0;
	for(i=0; i< (DeviceCell->length); i++)
	{
	  m= (i % 2); //which atom in 4 atom unit cell is i
	  
	  if(m==0)
	  {
	    if( (DeviceCell->siteinfo)[i][0] == 0)
	    {
	     l++; 
	    }
	  }
	}
	cellinfo->group_dim = l;
	//printf("#group %d\n", cellinfo->group_dim);

	cellinfo->group_sites = createIntArray(cellinfo->group_dim);
	l=0;
	for(i=0; i< (DeviceCell->length); i++)
	{
	  m= (i % 2); //which atom in 4 atom unit cell is i
	  
	  if(m==0)
	  {
	    if( (DeviceCell->siteinfo)[i][0] == 0)
	    {
	      (cellinfo->group_sites)[l] = i;
	      l++;
	    }
	  }
	}
	
      
      }
      
      
      
      
      
    }
    
    
    else if(config==2) //just the sites that actually connect to leads, and not the full conventional cells
		       //both leads here are in the starting cell, no group sites
    {
      j=0;
      cellinfo->num_leads=2;
      cellinfo->lead_dims=createIntArray(cellinfo->num_leads);
      
      if( (DeviceCell->geo) == 0)	//zigzag case
      {
	
	//LHS lead connections
	for(i=0; i< 2*(DeviceCell->length); i++)
	{
	  m= (i % 4); //which atom in 4 atom unit cell is i
	  
	  if(m==1 || m==2)
	  {
	 
	    if( (DeviceCell->siteinfo)[i][0] == 0)
	    {
		(cellinfo->cells_site_order)[j] = i;
		j++;
		
		(cellinfo->sites_by_cell)[i] = 0;
	    }
	
	  }
	}
	(cellinfo->lead_dims)[0] = j;
	
	
	//RHS lead connections
	for(i=0; i< 2*(DeviceCell->length); i++)
	{
	  m= (i % 4); //which atom in 4 atom unit cell is i
	  
	  if(m==0 || m==3)
	  {
	 
	    if( (DeviceCell->siteinfo)[alldim - 2*(DeviceCell->length) + i][0] == 0)
	    {
		(cellinfo->cells_site_order)[j] = alldim - 2*(DeviceCell->length) + i;
		j++;
		
		(cellinfo->sites_by_cell)[alldim - 2*(DeviceCell->length) + i] = 0;
	    }
	
	  }
	}
	cellinfo->cell1dim = j;  
	(cellinfo->lead_dims)[1] = j- (cellinfo->lead_dims)[0];
	
	//printf("#cell1 %d\n", cellinfo->cell1dim);
	
	//group sites
	
	cellinfo->group_dim = 0;
	
	cellinfo->lead_sites = createIntArray( (cellinfo->lead_dims)[0] + (cellinfo->lead_dims)[1]);
	
	for(l=0; l<(cellinfo->lead_dims)[0] + (cellinfo->lead_dims)[1]; l++)
	{
	  cellinfo->lead_sites[l] = (cellinfo->cells_site_order)[l];
	}
	
	
	
	//printf("#group %d\n", cellinfo->group_dim);

		
      
      }
      
      
      if( (DeviceCell->geo) == 1)	//armchair case
      {
	
	//LHS connections
	for(i=0; i< (DeviceCell->length); i++)
	{
	  m= (i % 2); //which atom in 2 atom unit cell is i
	  
	  if(m==0 )
	  {
	 
	    if( (DeviceCell->siteinfo)[i][0] == 0)
	    {
		(cellinfo->cells_site_order)[j] = i;
		j++;
		
		(cellinfo->sites_by_cell)[i] = 0;
	    }
	
	  }
	}
	(cellinfo->lead_dims)[0] = j;
	
	//RHS connections
	for(i=0; i< (DeviceCell->length); i++)
	{
	  m= (i % 2); //which atom in 2 atom unit cell is i
	  
	  if(m==0 )
	  {
	 
	    if( (DeviceCell->siteinfo)[alldim - (DeviceCell->length) + i][0] == 0)
	    {
		(cellinfo->cells_site_order)[j] = alldim -(DeviceCell->length) + i;
		j++;
		
		(cellinfo->sites_by_cell)[alldim - (DeviceCell->length) + i] = 0;
	    }
	
	  }
	}
	cellinfo->cell1dim = j;  
	(cellinfo->lead_dims)[1] = j -(cellinfo->lead_dims)[0];
	//printf("#cell1 %d\n", cellinfo->cell1dim);
		 // printf("lead dims %d:	%d\n", (cellinfo->lead_dims)[0],  (cellinfo->lead_dims)[1]);

	//group sites
	
	cellinfo->group_dim = 0;
	//printf("#group %d\n", cellinfo->group_dim);

	cellinfo->lead_sites = createIntArray( (cellinfo->lead_dims)[0] + (cellinfo->lead_dims)[1]);
	
	for(l=0; l<(cellinfo->lead_dims)[0] + (cellinfo->lead_dims)[1]; l++)
	{
	  cellinfo->lead_sites[l] = (cellinfo->cells_site_order)[l];
	  //printf("leadsite %d:	%d\n", l,  cellinfo->lead_sites[l]);
	}
	
	
      
      }
      
      
      
      
      
    }
    
    
    
    
    
//     else if(config==1) //test run with both sides
//     {
//       j=0;
//       
//       for(i=0; i< 2*(DeviceCell->length); i++)
//       {
// 	if( (DeviceCell->siteinfo)[i][0] == 0)
// 	{
// 	    (cellinfo->cells_site_order)[j] = i;
// 	    j++;
// 	    
// 	    (cellinfo->sites_by_cell)[i] = 0;
// 	}
// 	
//       }
//       
//       for(i=0; i< 2*(DeviceCell->length); i++)
//       {
// 	if( (DeviceCell->siteinfo)[alldim - 2*(DeviceCell->length) + i][0] == 0)
// 	{
// 	    (cellinfo->cells_site_order)[j] = alldim - 2*(DeviceCell->length) + i;
// 	    j++;
// 	    
// 	    (cellinfo->sites_by_cell)[alldim - 2*(DeviceCell->length) + i] = 0;
// 	}
// 	
//       }
//       
//       cellinfo->cell1dim = j;
//     }
    
    else
      exit(1);
  
}

void cellSplitter (RectRedux *DeviceCell, cnxProfile *cnxp, cellDivision *cellinfo)
{
  int i, j, k, l, m;
  
  
  int dim = *(DeviceCell->Nrem);
  int alldim =2*(DeviceCell->length)*(DeviceCell->length2);
  //int *sites_left = createIntArray(dim);
  //int *sites_still_left = createIntArray(dim);
  int num_sites_left=0, num_sites_still_left=0;
  int checker, cellindex, last_cell_start, last_cell_size, cell_start, cell_size, sites_to_date;
  int added_with_group=0, group_added=0;
  
  num_sites_left=0;
  cellinfo->group_cell=-1;
  
  
  for(i=0; i<alldim; i++)
  {
    if( (cellinfo->sites_by_cell)[i] == -1)
    {
      //sites_left[num_sites_left] = i;
      num_sites_left++;
    }
  }
  
  cellindex=0; last_cell_start=0; last_cell_size=(cellinfo->cell1dim); sites_to_date=(cellinfo->cell1dim);
  
  while(num_sites_left>0)
  {
    cellindex++;
    //num_sites_still_left=0; 
    cell_start = last_cell_start+last_cell_size;
    cell_size = 0;
   //printf("cell %d - added sites\t", cellindex);
    
    group_added=0;
    
    //add neighbours of sites in previous cell which are not yet added
    for(i=last_cell_start; i<last_cell_start + last_cell_size; i++)
    {
      k=(cellinfo->cells_site_order)[i]; 
      for(j=0; j<(cnxp->site_cnxnum)[k]; j++)
      {
	l=(cnxp->site_cnx)[k][j];
	
	if ((cellinfo->sites_by_cell)[l] == -1)
	{
	  added_with_group=0;
	  
	  //check if site just added is in the SE group, and if so add all other sites in this group
	  //group addition seems to work fine!
	  if(group_added == 0 && (cellinfo->group_dim) > 0)
	  {
	    for(m=0; m<(cellinfo->group_dim); m++)
	    {
	      if(l==(cellinfo->group_sites[m]))
	      {
		added_with_group=1;
		group_added=1;
	      }
	    }
	    
	    if(group_added == 1)
	    {
	      cellinfo->group_cell = cellindex;
	      for(m=0; m<(cellinfo->group_dim); m++)
	      {
		(cellinfo->cells_site_order)[sites_to_date] = (cellinfo->group_sites[m]);
		sites_to_date++;
		num_sites_left--;
		(cellinfo->sites_by_cell)[(cellinfo->group_sites[m])] = cellindex;
		cell_size++;
	      }
	    }
	  }
	 
	  //else just add this site
	  if(added_with_group==0)
	  {
	    (cellinfo->cells_site_order)[sites_to_date] = l;
	    sites_to_date++;
	    num_sites_left--;
	    (cellinfo->sites_by_cell)[l] = cellindex;
	    cell_size++;
	    //printf("%d\t", l);
	  }
	 
	  
	}	
      }
      
    }
    
    //printf("\n");
    last_cell_size=cell_size;
    last_cell_start=cell_start;
        
  }
    
  (cellinfo->num_cells)=cellindex+1;
  cellinfo->starting_index = createIntArray((cellinfo->num_cells));
  cellinfo->cell_dims = createIntArray((cellinfo->num_cells));
  
  l=0;
  for(i=0; i< (cellinfo->num_cells); i++)
  {
    (cellinfo->starting_index)[i] = l;
    k=0;
    for(j=0; j<alldim; j++)
    {
      if((cellinfo->sites_by_cell)[j] == i)
      {
	k++;
      }
    }
    (cellinfo->cell_dims)[i] = k;
    l += k;
  }
      
//       for(i=0;i<(cellinfo->num_cells);i++)
//       {
// 	printf("cell %d : ", i);
// 	for(j=0; j<(cellinfo->cell_dims)[i]; j++)
// 	{
// 	  printf("%d\t", (cellinfo->cells_site_order)[(cellinfo->starting_index)[i] + j]);
// 	}
// 	printf("\n");
//       }

//     for(i=0;i<(cellinfo->num_cells);i++)
//     {
//       for(j=0; j<(cellinfo->cell_dims)[i]; j++)
//       {
// 	printf("%lf	%lf\n", (DeviceCell->pos)[(cellinfo->cells_site_order)[(cellinfo->starting_index)[i] + j]][0], (DeviceCell->pos)[(cellinfo->cells_site_order)[(cellinfo->starting_index)[i] + j]][1]);
//       }
//       printf("\n");
//     }
    
  
}




//This function simply generates the connections present in a pristine nearest neighbour AC or ZZ cell
//Require a DeviceCell using the "geo" notation
//says whether i and j are connected or not
int zzacnn (RectRedux *DeviceCell, void *rule_params, int a, int b)
{
    int i, j, k, l, ans;  //set ans =0 if there is a connection
    ans=1;
    
    
      //zigzag case - currently failing for odd-index ZGNRs
      if ((DeviceCell->geo) == 0)
      {
	      i= (a % (2 * (DeviceCell->length))); 	//which atom in chain is a?
	      j= (i % 4);				//and which atom in 4 atom unit cell is it?
	      
	      if(j==0)
	      {
		if(b==(a-1) && (i!=0))
		{
		  ans=0;
		}
		
		if(b==(a+1))
		{
		  ans=0;
		}
		
		if(b==(a+1 + 2 *(DeviceCell->length)))
		{
		  ans=0;
		}
	      }
	      
	      if(j==1)
	      {
		if(b==(a-1))
		{
		  ans=0;
		}
		
		if(b==(a+1) && (i!=2 *(DeviceCell->length) -1))
		{
		  ans=0;
		}
		
		if(b==(a-1 - 2 *(DeviceCell->length)))
		{
		  ans=0;
		}
	      }
	      
	      if(j==2)
	      {
		if(b==(a-1))
		{
		  ans=0;
		}
		
		if(b==(a+1))
		{
		  ans=0;
		}
		
		if(b==(a+1 - 2 *(DeviceCell->length)))
		{
		  ans=0;
		}
	      }
	      
	      if(j==3)
	      {
		if(b==(a-1) )
		{
		  ans=0;
		}
		
		if(b==(a+1)&& (i!=2 *(DeviceCell->length) -1))
		{
		  ans=0;
		}
		
		if(b==(a-1 + 2 *(DeviceCell->length)))
		{
		  ans=0;
		}
	      }
      }
      
      //armchair case
      if((DeviceCell->geo) == 1)
      {
	  i= (a % (2 * (DeviceCell->length))); //which atom in unit cell
	  k= (a % ( (DeviceCell->length)));	//which atom in chain
	  j= (a / (2 * (DeviceCell->length))); //which unit cell
	  l= (i / (DeviceCell->length)); 	//which chain in unit cell
	  
	  
	  if(l==0)
	  {
	      if( (k%2) == 0)
	      {
		if(b==a-(DeviceCell->length))
		{
		  ans=0;
		}
		if(b==a-1)
		{
		  if(k!=0)
		  {
		    ans=0;
		  }
		}
		if(b==a+1)
		{
		  if(k!= (DeviceCell->length)-1 )
		  {
		    ans=0;
		  }
		}
	      }
	      
	      if( (k%2) == 1)
	      {
		if(b==a+(DeviceCell->length))
		{
		  ans=0;
		}
		if(b==a-1)
		{
		  if(k!=0)
		  {
		    ans=0;
		  }
		}
		if(b==a+1)
		{
		  if(k!= (DeviceCell->length)-1 )
		  {
		    ans=0;
		  }
		}
	      }
	  }
	  
	  if(l==1)
	  {
	      if( (k%2) == 0)
	      {
		if(b==a+(DeviceCell->length))
		{
		  ans=0;
		}
		if(b==a-1)
		{
		  if(k!=0)
		  {
		    ans=0;
		  }
		}
		if(b==a+1)
		{
		  if(k!= (DeviceCell->length)-1 )
		  {
		    ans=0;
		  }
		}
	      }
	      
	      if( (k%2) == 1)
	      {
		if(b==a-(DeviceCell->length))
		{
		  ans=0;
		}
		if(b==a-1)
		{
		  if(k!=0)
		  {
		    ans=0;
		  }
		}
		if(b==a+1)
		{
		  if(k!= (DeviceCell->length)-1 )
		  {
		    ans=0;
		  }
		}
	      }
	  }
	  
	
      }

    
    return ans;
}



//This function simply generates the connections present in a pristine nearest neighbour AC or ZZ cell
//Require a DeviceCell using the "geo" notation
//says whether i and j are connected or not
//this version includes ky periodicity - but only through connectivity - does not say what Bloch phase is
int zzacnnk (RectRedux *DeviceCell, void *rule_params, int a, int b)
{
    int i, j, k, l, ans;  //set ans =0 if there is a connection
    ans=1;
    
    
      //zigzag case 
      if ((DeviceCell->geo) == 0)
      {
	      i= (a % (2 * (DeviceCell->length))); 	//which atom in chain is a?
	      j= (i % 4);				//and which atom in 4 atom unit cell is it?
	      
	      if(j==0)
	      {
		if(b==(a-1) && (i!=0))
		{
		  ans=0;
		}
		
		if(b==(a+1))
		{
		  ans=0;
		}
		
		if(b==(a+1 + 2 *(DeviceCell->length)))
		{
		  ans=0;
		}
	      }
	      
	      if(j==1)
	      {
		if(b==(a-1))
		{
		  ans=0;
		}
		
		if(b==(a+1) && (i!=2 *(DeviceCell->length) -1))
		{
		  ans=0;
		}
		
		if(b==(a-1 - 2 *(DeviceCell->length)))
		{
		  ans=0;
		}
	      }
	      
	      if(j==2)
	      {
		if(b==(a-1))
		{
		  ans=0;
		}
		
		if(b==(a+1))
		{
		  ans=0;
		}
		
		if(b==(a+1 - 2 *(DeviceCell->length)))
		{
		  ans=0;
		}
	      }
	      
	      if(j==3)
	      {
		if(b==(a-1) )
		{
		  ans=0;
		}
		
		if(b==(a+1)&& (i!=2 *(DeviceCell->length) -1))
		{
		  ans=0;
		}
		
		if(b==(a-1 + 2 *(DeviceCell->length)))
		{
		  ans=0;
		}
	      }
	      
	      //the periodicity ones!
	      if(i== (2 *(DeviceCell->length) -1))
	      {
		if(b== (a - 2 *(DeviceCell->length) + 1))
		{
		  ans=0;
		}
	      }
	      
	      if(i==0)
	      {
		if(b== a+ 2 *(DeviceCell->length) -1)
		{
		  ans=0;
		}
	      }
      }
      
      //armchair case
      if((DeviceCell->geo) == 1)
      {
	  i= (a % (2 * (DeviceCell->length))); //which atom in unit cell
	  k= (a % ( (DeviceCell->length)));	//which atom in chain
	  j= (a / (2 * (DeviceCell->length))); //which unit cell
	  l= (i / (DeviceCell->length)); 	//which chain in unit cell
	  
	  
	  if(l==0)
	  {
	      if( (k%2) == 0)
	      {
		if(b==a-(DeviceCell->length))
		{
		  ans=0;
		}
		if(b==a-1)
		{
		  if(k!=0)
		  {
		    ans=0;
		  }
		}
		if(b==a+1)
		{
		  if(k!= (DeviceCell->length)-1 )
		  {
		    ans=0;
		  }
		}
	      }
	      
	      if( (k%2) == 1)
	      {
		if(b==a+(DeviceCell->length))
		{
		  ans=0;
		}
		if(b==a-1)
		{
		  if(k!=0)
		  {
		    ans=0;
		  }
		}
		if(b==a+1)
		{
		  if(k!= (DeviceCell->length)-1 )
		  {
		    ans=0;
		  }
		}
	      }
	  }
	  
	  if(l==1)
	  {
	      if( (k%2) == 0)
	      {
		if(b==a+(DeviceCell->length))
		{
		  ans=0;
		}
		if(b==a-1)
		{
		  if(k!=0)
		  {
		    ans=0;
		  }
		}
		if(b==a+1)
		{
		  if(k!= (DeviceCell->length)-1 )
		  {
		    ans=0;
		  }
		}
	      }
	      
	      if( (k%2) == 1)
	      {
		if(b==a-(DeviceCell->length))
		{
		  ans=0;
		}
		if(b==a-1)
		{
		  if(k!=0)
		  {
		    ans=0;
		  }
		}
		if(b==a+1)
		{
		  if(k!= (DeviceCell->length)-1 )
		  {
		    ans=0;
		  }
		}
	      }
	  }
	  
	     //the periodicity ones!
	  if(k==(DeviceCell->length)-1)
	  {
	    if(b== (a-(DeviceCell->length)+1))
	    {
	      ans=0;
	    }
	    
	  }
	  if(k==0)
	  {
	    if(b== (a+(DeviceCell->length)-1))
	    {
	      ans=0;
	    }
	    
	  }
      
      }
      
 

    
    return ans;
}

