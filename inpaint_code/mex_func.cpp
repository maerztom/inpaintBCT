/* mex_func.cpp  --- inpaintBCT 
 * Copyright (C) 2013 Thomas März (maerz@maths.ox.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "global_def.h"

#define NO_ERR 0
#define ERR_NARGIN 1
#define ERR_DOUBLE 2
#define ERR_STRING 3
#define ERR_DIM 4

#define ERR_PARAM_VAL_N 5
#define ERR_PARAM_VAL_C 6
#define ERR_PARAM_DIM_N 7
#define ERR_PARAM_DIM_C 8
#define ERR_NO_MASK_ORDER 9
#define ERR_COEFF 10
#define ERR_ORDER 11
#define ERR_ARG_MISSING 12 
#define ERR_UNKNOWN_ID 13
#define ERR_MASK_DIM 14
#define ERR_EMPTY_MASK 15

#define TYPE_N 0
#define TYPE_C 1
#define TYPE_D 2

void *AllocMem(size_t n)
{
	void *p;
	p = mxMalloc( n );
	mexMakeMemoryPersistent(p);
	return p;
}

void FreeMem(void* p)
{
	mxFree(p);
}

void ErrorMessage(int type)
{
    switch( type )
    {
        case ERR_NARGIN:
            mexPrintf("Error: The number of input arguments must be odd and greater than 2! \n");
            break;
            
        case ERR_DOUBLE:
            mexPrintf("Error: Data type of all numeric arguments must be double \n");
            break;
            
        case ERR_DIM:
            mexPrintf("Error: Image must be n x m matrix or n x m x k multi-matricx for k-channel images \n");
            break;
            
        case ERR_MASK_DIM:
            mexPrintf("Error: Mask <-> Image dimesion mismatch \n");
            break;
            
        case ERR_EMPTY_MASK:
            mexPrintf("Error: Empty mask, nothing to do \n");
            break;
            
        case ERR_ARG_MISSING:
            mexPrintf("Error: Argument missing \n");
            break;
            
        case ERR_UNKNOWN_ID:
            mexPrintf("Error: Unknown specifier \n");
            break;
            
        case ERR_STRING:
            mexPrintf("Error: All argument specifiers must be strings \n");
            break;
            
        case ERR_PARAM_VAL_N:
            mexPrintf("Epsilon must be greater than or equal 1. \n");
			break;
            
        case ERR_PARAM_VAL_C:
            mexPrintf("Error in the parameter values: \n");
            mexPrintf("Epsilon must be greater than or equal 1. \n");
            mexPrintf("The parameters kappa and sigma must be greater than or equal zero. \n");
			mexPrintf("Rho must be greater than zero \n");
			mexPrintf("Thresh must be greater than or equal zero \n");
            mexPrintf("Quant must be greater than zero \n");
			break;
            
        case ERR_PARAM_DIM_N:
            mexPrintf("Error: Parameter must be a scalar \n");
            break;
            
        case ERR_PARAM_DIM_C:
            mexPrintf("Error: Parameter list must be a vector \n");
            break;
            
        case ERR_COEFF:
            mexPrintf("Error in the channel coefficients: \n");
            mexPrintf("The channel coefficients must be a vector and \n");
            mexPrintf("the dimension of the coefficients vector must match the number of channels. \n");
            mexPrintf("Any coefficient must be greater than or equal zero and the sum of the coefficients must be greater than zero. \n");
            break;
            
        case ERR_NO_MASK_ORDER:
            mexPrintf("Error: Mask has to be specified, there is no default for this parameter. \n");
            break;
            
            
        case ERR_ORDER:
            mexPrintf("Error: The combined mask and order data M must be a 3 x N matrix, \n");
            mexPrintf("       where M(:,k) = [i(k) ; j(k) ; T(i,j)] and N is the number of points to be inpainted. \n");
            break;
    }
}

void SetDefaults(Data *data)
{
    data->rows = 1;
    data->cols = 1;
    data->channels = 1;
    data->size = 1;
    data->Image = NULL;
    data->MImage = NULL;

    

    // default parameters
    data->epsilon = 5;
    data->radius = 5;
    data->kappa = 25;
    data->sigma = 1.414213562373095; // sqrt(2.0);
    data->rho = 5;
    data->thresh = 0;
    data->delta_quant4 = 1; // default quantization range is [0,255]

    data->convex = NULL;


    data->ordered_points = NULL;
    data->nof_points2inpaint = 0;
    data->heap = NULL;
    data->Tfield = NULL;
    data->Domain = NULL;
    data->MDomain = NULL;
    data->GivenGuidanceT = NULL;
    
    data->lenSK1 = 0;
    data->lenSK2 = 0;
    data->SKernel1 = NULL;
    data->SKernel2 = NULL;
    data->Ihelp = NULL;
    data->Shelp = NULL;
    
    data->ordergiven = 0;
    data->guidance = 1;
    
    data->inpaint_undefined = 0;
}
    
void ClearMemory(Data *data)
{
    if( data->Domain != NULL )
    {
        FreeMem( data->Domain );
        data->Domain = NULL;
    }
    
    if( data->MDomain != NULL )
    {
        FreeMem( data->MDomain );
        data->MDomain = NULL;
    }
    
    if( data->Tfield != NULL )
    {
        FreeMem( data->Tfield );
        data->Tfield = NULL;
    }
    
    if( data->heap != NULL )
    {
        FreeMem( data->heap );
        data->heap = NULL;
    }
    
    if( data->Ihelp != NULL )
    {
        FreeMem( data->Ihelp );
        data->Ihelp = NULL;
    }
    
    if( data->convex != NULL )
    {
        FreeMem( data->convex );
        data->convex = NULL;
    }
    
    if( data->SKernel1 != NULL )
    {
        FreeMem( data->SKernel1 );
        data->SKernel1 = NULL;
    }
    
    if( data->SKernel2 != NULL )
    {
        FreeMem( data->SKernel2 );
        data->SKernel2 = NULL;
    }
    
    if( data->Shelp != NULL )
    {
        FreeMem( data->Shelp );
        data->Shelp = NULL;
    }
    
    if( data->ordered_points != NULL )
    {
        FreeMem( data->ordered_points );
        data->ordered_points = NULL;
    }     
}

int GetMask( const mxArray *arg, Data *data)
{
    int err = NO_ERR;
    int nof_dim;
    double *parg;
    int i,j,c;
    int index;
    int not_equal;
    
    
    if( !mxIsDouble(arg) ) 
    {
        err = ERR_DOUBLE;
        return err;
    }
    
    
    if( (mxGetM( arg ) != data->rows) || (mxGetN( arg ) != data->cols) )
    {
        err = ERR_MASK_DIM; 
        return err;
    }
    
    data->nof_points2inpaint = 0;
    parg = mxGetPr( arg );
    
    for( i = 0; i < data->rows ; i++)
    {
        for( j = 0; j < data->cols ; j++)
        {
            index = j * data->rows + i;

            if( parg[index] != 0 )
            {
                data->nof_points2inpaint = data->nof_points2inpaint + 1;
                data->Domain[index] = 0; //INSIDE
                data->MDomain[index] = 0;

                for(c = 0; c < data->channels ; c++)
                {
                    data->Image[index + c * data->size] = 0;
                    data->MImage[index + c * data->size] = 0;
                }
            }
            else
            {
                data->Domain[index] = 1; // OUTSIDE
                data->MDomain[index] = 1;
            }

        }
    }
    
    if( data->nof_points2inpaint == 0 )
        err = ERR_EMPTY_MASK;
    
    return err;
}

int GetOrder( const mxArray *arg, Data *data)
{
    int err = NO_ERR;
    int rows;
    int size;
    double *parg;
    int i,j,c;
    int index;
    
    if( !mxIsDouble(arg) ) 
    {
        err = ERR_DOUBLE;
        return err;
    }
    
    rows = mxGetM( arg );
    if( rows != 3)
    {
        err = ERR_ORDER;
        return err;
    }

    // transfer info about ordered points
    data->nof_points2inpaint = mxGetN( arg );
    
    size = 3 * data->nof_points2inpaint;     
    parg = mxGetPr( arg );

    for( i=0 ; i < size ; i++ )
        if( (i+1) % 3 != 0 )
            data->ordered_points[i] = parg[i] - 1;
        else
            data->ordered_points[i] = parg[i];

    // transfer Tfield
    for( j = 0; j < data->cols ; j++)
    {
        for( i = 0; i < data->rows ; i++)
        {
            index = j*data->rows + i;
            
            // OUTSIDE
            data->Domain[index] = 1; 
            data->MDomain[index] = 1;

            data->Tfield[index].T = -1;
            data->Tfield[index].flag = KNOWN;
            data->Tfield[index].hpos = -1;
            data->Tfield[index].i = i;
            data->Tfield[index].j = j;
        }
    }
    
    for( i=0 ; i < size ; i=i+3 )
    {
        index = data->ordered_points[i+1] * data->rows + data->ordered_points[i];

        // INSIDE
        data->Domain[index] = 0; 
        data->MDomain[index] = 0;

        for(c = 0; c < data->channels ; c++)
        {
            data->Image[index + c * data->size] = 0;
            data->MImage[index + c * data->size] = 0;
        }

        data->Tfield[index].T = data->ordered_points[i+2];
        data->Tfield[index].flag = TO_INPAINT;
        data->Tfield[index].hpos = -1;
        data->Tfield[index].i = data->ordered_points[i];
        data->Tfield[index].j = data->ordered_points[i+1];
    }
    
    return err;
}

void SetKernels(Data *data)
{
    int i;
	int s;
	int r;
    
	s = max( round(2 * data->sigma) , 1 );
	r = max( round(2 * data->rho) , 1 ); 
	data->lenSK1 = 2*s +1;
	data->lenSK2 = 2*r +1;

    
    if( data->sigma > 0 )
    {
        data->SKernel1 = (double *)AllocMem(sizeof(double) * data->lenSK1);
        for( i=0 ; i < data->lenSK1 ; i++)
            data->SKernel1[i] = exp( -((i-s)*(i-s))/(2* data->sigma * data->sigma) );
        
        data->Shelp = (double *) AllocMem(sizeof(double) * data->lenSK1);
    }
    
    data->SKernel2 = (double *)AllocMem(sizeof(double) * data->lenSK2);
    for( i=0 ; i < data->lenSK2 ; i++)
        data->SKernel2[i] = exp( -((i-r)*(i-r))/(2* data->rho * data->rho) );
    
}

int GetParam( const mxArray *arg, Data *data , int type)
{
    int err = NO_ERR;
    int nof_param;
	double *paramlist;
    
    if( !mxIsDouble(arg) ) 
    {
        err = ERR_DOUBLE;
        return err;
    }
    
    nof_param = max(mxGetM(arg),mxGetN(arg));
    paramlist = mxGetPr(arg);
        
    if( type == TYPE_N)
    {
        if( nof_param != 1 )
        {
            err = ERR_PARAM_DIM_N;
            return err;
        }

        data->epsilon = paramlist[0];
        data->radius = (int) (data->epsilon + 0.5);
        
        if( data->epsilon < 1 )
        {
            err = ERR_PARAM_VAL_N;
            return err;
        }
    }
    else if( type == TYPE_C)
    {
        if( min(mxGetM(arg),mxGetN(arg)) != 1 )
        {
            err = ERR_PARAM_DIM_C;
            return err;
        }

        if(nof_param >= 1)
        {
            data->epsilon = paramlist[0];
            data->radius = (int) (data->epsilon + 0.5);
        }
        if(nof_param >= 2)
            data->kappa = paramlist[1];
        if(nof_param >= 3)
            data->sigma = paramlist[2];
        if(nof_param >= 4)
            data->rho = paramlist[3];
        if(nof_param >= 5)
            data->thresh = paramlist[4];
        if(nof_param >= 6)
        {
            //data->delta_quant4 = paramlist[5];
            data->delta_quant4 = paramlist[5]/255;
            data->delta_quant4 = data->delta_quant4 * data->delta_quant4;
            data->delta_quant4 = data->delta_quant4 * data->delta_quant4;
        }
        
        // check parameters 
        if( (data->epsilon < 1 ) || ( data->kappa < 0 ) || ( data->sigma < 0) || ( data->rho <= 0 ) || ( data->thresh < 0) || (data->delta_quant4 == 0) )
        {
            err = ERR_PARAM_VAL_C;
            return err;
        }
        
        if(nof_param >= 7)
        {
            double sum;
            int c;
            
            if(nof_param != (6 + data->channels))
            {
                err = ERR_COEFF;
                return err;
            }
            
            data->convex = (double *)AllocMem(sizeof(double) * data->channels);
            
            for( c=0 ; c < data->channels ; c++)
            {
                if( paramlist[6+c] < 0 )
                {
                    err = ERR_COEFF;
                    return err;
                }
                sum = sum + paramlist[6+c];
            }
            
            for( c=0 ; c < data->channels ; c++)
                data->convex[c] = paramlist[6+c]/sum;
        }
    }
    else
    {
        data->GivenGuidanceT = paramlist;
    }
    return err;
}






// debugging procs
// void CopyMask( Data *data )
// {
//     int i,j,index;
//     
//     for( i = 0; i < data->rows ; i++)
//     {
//         for( j = 0; j < data->cols ; j++)
//         {
//             index = j * data->rows + i;
//             data->Image[index] = data->MDomain[index];
//         }
//     }
// }
// 
// void CopyParam( Data *data )
// {
//     int i;
//     data->Image[0] = data->radius;
//     data->Image[1] = data->epsilon;
//     data->Image[2] = data->kappa;
//     data->Image[3] = data->sigma;
//     data->Image[4] = data->rho;
//     data->Image[5] = data->thresh;
//     data->Image[6] = data->delta_quant4;
//     
//     if( data->convex != NULL )
//         for( i=0 ; i < data->channels ; i++ )
//             data->Image[7+i] = data->convex[i];
//     
// }
// 
// void CopyDirfield( Data *data )
// {
//     int i,j,c,index;
//     
//     if( data->DirField == NULL )
//         return;
//     
//     
//     for( i = 0; i < data->rows ; i++)
//     {
//         for( j = 0; j < data->cols ; j++)
//         {
//             index = j * data->rows + i;
//             
//             for( c=0 ; c < data->channels ; c++ )
//                 data->Image[index + c * data->size] = data->DirField[index + c * data->size];
//         }
//     }
// }
// 
// void CopyMImage( Data *data )
// {
//     int i,j,c,index;
//     
//     for( i = 0; i < data->rows ; i++)
//     {
//         for( j = 0; j < data->cols ; j++)
//         {
//             index = j * data->rows + i;
//             
//             if( data->MDomain[index] > 1e-3 )
//                 for( c = 0 ; c < data->channels ; c++ )
//                     data->Image[index + c * data->size] = data->MImage[index + c * data->size]/data->MDomain[index];
//         }
//     }
// }
// 
// void CopyTfield( Data *data )
// {
//     int i,j,index;
//     
//     for( i = 0; i < data->rows ; i++)
//     {
//         for( j = 0; j < data->cols ; j++)
//         {
//             index = j * data->rows + i;
//             data->Image[index] = data->Tfield[index].T;
//         }
//     }
// }
// end debuging procs



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Data data;
	mxArray *Im;
	int nof_dim;
    const int *elem_per_dim;
    int nargout;
    int i;
	
    
    // set default values
    SetDefaults(&data);
    
	// check number of input arguments
	if( (nrhs % 2 == 0) || (nrhs < 3 ) )
	{
		ErrorMessage(ERR_NARGIN);
        ClearMemory(&data);
		return;
	}
    
	// initialize that part of the data depending only on input image
    {
        if( !mxIsDouble(prhs[0]) ) // check if image is given as double
        {
            ErrorMessage(ERR_DOUBLE);
            ClearMemory(&data);
            return;
        }
        
        // get image dimensions
        nof_dim = mxGetNumberOfDimensions(prhs[0]);
        elem_per_dim = mxGetDimensions(prhs[0]);
        
        if( nof_dim > 3 )
        {
            ErrorMessage(ERR_DIM);
            ClearMemory(&data);
            return;
        }
        
        plhs[0] = mxDuplicateArray(prhs[0]); // output image
        Im = mxDuplicateArray(prhs[0]);      // smoothed image

        data.rows = elem_per_dim[0];
        data.cols = elem_per_dim[1];
        data.size = data.rows * data.cols;
        data.channels = 1;
        if( nof_dim == 3 )
            data.channels = elem_per_dim[2];
        
        data.Image = mxGetPr(plhs[0]);
        data.MImage = mxGetPr(Im); 
        
        data.Ihelp = (double *)AllocMem(sizeof(double) * data.channels);
        data.Tfield = (hItem *) AllocMem(sizeof(hItem) * data.size);
        data.Domain = (double *) AllocMem(sizeof(double) * data.size);
        data.MDomain = (double *) AllocMem(sizeof(double) * data.size);
        data.heap = (hItem **) AllocMem(sizeof(hItem *) * data.size);
        data.ordered_points = (double *) AllocMem(sizeof(double) * data.size *3);
        
        // number of output arguments
        if( nlhs <= 2 )
            nargout = nlhs;
        else
            nargout = 2;
    }
     
    
    // get and check arguments
    {
        char guidance_done = 0;
        char order_done = 0;
        char argname[100];
        int  len;
        int  err;
        
        for( i = 1 ; i < nrhs ; i = i+2 )
        {
            if( !mxIsChar( prhs[i] ) )
            {
                ErrorMessage(ERR_STRING);
                ClearMemory(&data);
                return;
            }
            else
            {
                len = mxGetN( prhs[i] );
                mxGetString( prhs[i], argname, len+1);
                
                //mexPrintf(" %s \n", argname);
                
                switch( argname[0] )
                {
                    case 'o': // orderD , orderT
                        {
                            if( !order_done )
                            {
                                if( !mxIsEmpty(prhs[i+1]) )
                                {
                                    if( argname[len-1] == 'D' )
                                    {
                                        err = GetMask( prhs[i+1], &data);
                                        data.ordergiven = 0;
                                    }
                                    else if( argname[len-1] == 'T' )
                                    {
                                        err = GetOrder( prhs[i+1], &data);
                                        data.ordergiven =1;
                                    }
                                    else
                                        err = ERR_UNKNOWN_ID;
                                }
                                else
                                    err = ERR_ARG_MISSING; 

                                order_done = 1;
                            }
                            break;
                        }
                                            
                    case 'g': // guidanceN , guidanceC , guidanceD
                        {
                            if( !guidance_done )
                            {
                                if( !mxIsEmpty(prhs[i+1]) )
                                {
                                    if( argname[len-1] == 'N' )
                                    {
                                        err = GetParam( prhs[i+1], &data, TYPE_N);
                                        data.guidance = 0;
                                    }
                                    else if( argname[len-1] == 'C' )
                                    {
                                        err = GetParam( prhs[i+1], &data, TYPE_C);
                                        data.guidance = 1;
                                    }
                                    else if( argname[len-1] == 'D' )
                                    {
                                        err = GetParam( prhs[i+1], &data, TYPE_D);
                                        data.guidance = 2;
                                    }
                                    else
                                        err = ERR_UNKNOWN_ID;
                                }
                                else
                                    err = ERR_ARG_MISSING;
                                
                                // guidance_done = 1;
                                if(data.GivenGuidanceT != NULL)
                                    data.guidance = 2;
                            }
                            break;
                        }
                        
                    default:
                        err = ERR_UNKNOWN_ID;
                        break;
                }
                    
                if( err )
                {
                    ErrorMessage(err);
                    ClearMemory(&data);
                    return;
                }
            }
        }
        
        if( !order_done ) // at least mask must be specified
        {
            ErrorMessage(ERR_NO_MASK_ORDER);
            ClearMemory(&data);
            return;
        }
    }
	
    
    if( data.guidance == 1)
        SetKernels(&data);
    
	InpaintImage(&data); 
    
    if(data.inpaint_undefined == 1)
    {
        mexPrintf("\n\n");
        mexPrintf("Error:\n");
        mexPrintf("Some inpainted image values are undefined !\n");
        mexPrintf("This happens if the order is not well-defined. \n");
        mexPrintf("You can find out the undefined pixels by: > ind = find( isnan(result) ) \n\n\n");
    }
        
        
    if( nargout == 2 )
    {
        double *p;
        int size;
        
        size = 3 * data.nof_points2inpaint;
        plhs[1] = mxCreateDoubleMatrix(3, data.nof_points2inpaint , mxREAL);
        p = mxGetPr(plhs[1]);
        
        for( i=0 ; i < size ; i++ )
            if( (i+1) % 3 != 0 )
                p[i] = data.ordered_points[i] + 1;
            else
                p[i] = data.ordered_points[i];
    }
    
    //CopyMask(&data);
    //CopyParam( &data );
    //CopyDirfield( &data );
    //CopyMImage( &data );
    //CopyTfield(&data);
    
    ClearMemory(&data);
	return;
}

