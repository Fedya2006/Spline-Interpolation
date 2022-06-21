﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace ClassLibrary1
{
    //Класс для расчета значений функции в узлах
    public class Calculation
    {
        public static float ToSingle(string s)
        {
            return Convert.ToSingle(s);
        }

        public double[] Calc(int n, double[] x, Spf func)
        {
            double[] y = new double[n];
            switch (func)
            {
                case Spf.Cubic_polynomial:
                    for(int i = 0; i < n; i++)
                    {
                        y[i] = 0.1 * x[i] * x[i] * x[i] + 2 * x[i] * x[i] + 9 * x[i] + 1;
                       
                    }
                    return y;

                case Spf.Sin:
                    for (int i = 0; i < n; i++)
                    {
                        y[i] = Math.Sin(x[i]) + 6;
                    }
                    return y;

                case Spf.Random_func:
                    Random rand = new Random();
                    for (int i = 0; i < n; i++)
                    {
                        y[i] = rand.NextDouble() * 1000;
                    }
                    return y;
                default:
                    return new double[n];

            }
        }
    }
}
