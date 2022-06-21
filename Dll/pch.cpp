// pch.cpp: файл исходного кода, соответствующий предварительно скомпилированному заголовочному файлу

#include "pch.h"

// При использовании предварительно скомпилированных заголовочных файлов необходим следующий файл исходного кода для выполнения сборки.
extern "C" _declspec(dllexport)
void Spline_derivatives(double* x, double* y, int non_uniform_points_count, int uniform_points_count,  double left_derivative, double right_derivative, double* res, double& err)
{
    try
    {
        MKL_INT dorder[1] = {1};
        double* scoeff = new double[(non_uniform_points_count - 1) * DF_PP_CUBIC];
        DFTaskPtr task = new DFTaskPtr();

        int status = dfdNewTask1D(&task, non_uniform_points_count, x, DF_NON_UNIFORM_PARTITION, 1, y, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            err = 1;
            return;
        }
        double derivatives[] = { left_derivative , right_derivative };
        status = dfdEditPPSpline1D(task, DF_PP_CUBIC, DF_PP_NATURAL, DF_BC_2ND_LEFT_DER + DF_BC_2ND_RIGHT_DER, derivatives, DF_NO_IC, NULL, scoeff, DF_NO_HINT);
        if (status != DF_STATUS_OK)
        {
            err = 2;
            return;
        }
        status = dfdConstruct1D(task, DF_PP_SPLINE, DF_METHOD_STD);
        if (status != DF_STATUS_OK)
        {
            err = 3;
            return;
        }
        double range[] = { x[0], x[non_uniform_points_count - 1] };
        status = dfdInterpolate1D(task, DF_INTERP, DF_METHOD_PP, uniform_points_count, range, DF_UNIFORM_PARTITION, 1, dorder, NULL, res, DF_NO_HINT, NULL);
        if (status != DF_STATUS_OK)
        {
            err = 4;
            return;
        }
        status = dfDeleteTask(&task);
        if (status != DF_STATUS_OK)
        {
            err = 5;
            return;
        }
        
    }
    catch (int e)
    {
        err = e;
    }

}