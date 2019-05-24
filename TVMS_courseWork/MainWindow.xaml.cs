using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using MathNet;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Double;
using Microsoft.Win32;

namespace TVMS_courseWork
{
    /// <summary>
    /// Логика взаимодействия для MainWindow.xaml
    /// </summary>
    public partial class MainWindow : System.Windows.Window
    {
        double[][] columArray;
        double[][] columArraySource;

        double k = 43.2;
        int interval = 31;
        DenseMatrix koeffPair;
        DenseMatrix koeffPrivate;
        int colum;
        double t_tabl = 1.9640;
        double t_tabl2 = 1.15;
        double F_tabl = 1.755;
        double[] B;
        double incurracy;
        public MainWindow()
        {
            InitializeComponent();
            BtnDATA_Click();
            Normaliz();
            Pair();
            PrivateK();
            MultyplyK();
           // err();
            Regress();
        }
        #region Распределение
        private void Normaliz()
        {
            try
            {

                double[][] points = new double[colum][];
                double[][] xn = new double[colum][];
                double[][] m_e = new double[colum][];
                double[][] m_t = new double[colum][];
                double[] x2_arr = new double[colum];

                for (int i = 0; i < colum; i++)
                {
                    points[i] = point(columArray[i], interval);
                    double[] me_buf = m_e[i] = Get_m_e(columArray[i], points[i], interval);
                    xn[i] = Get_xn(columArray[i], points[i], interval);
                    double[] mt_buf = m_t[i] = theoretical_freq(columArray[i], m_e[i].Sum(), xn[i], interval);

                    for (int j = 0; j < mt_buf.Length; j++)
                    {
                        x2_arr[i] += Math.Pow(me_buf[j] - mt_buf[j], 2) / mt_buf[j];
                    }
                }
                string s = "\n", t = "\n", w = "\n\n";
                for (int i = 0; i < colum; i++)
                {
                    if (x2_arr[i] > k) { s += "Не нормальное\n\n"; t += string.Format("{0:F2}", x2_arr[i]) + "\n\n"; }
                    else { s += "Нормальное\n\n"; t += string.Format("{0:F2}", x2_arr[i]) + "\n\n"; }
                    w += Data.parametrs[i] + "\n\n";
                }

                stpDistr.Children.Add(new TextBlock { Text = s, HorizontalAlignment = HorizontalAlignment.Center });
                stpDistrX_2.Children.Add(new TextBlock { Text = t, HorizontalAlignment = HorizontalAlignment.Center });
                stpDistrParams.Children.Add(new TextBlock { Text = w });
            }
            catch (Exception ex) { MessageBox.Show(ex.Message, "Ошибка"); }
        }
        private double[] theoretical_freq(double[] arr, double sum, double[] xn, int CountInterval)
        {
            double X_aver = DiscriptiveStatistics.Average(arr);
            double disp = Math.Sqrt(DiscriptiveStatistics.Dispersion(arr));
            double[] m_t = new double[CountInterval];
            double h = (DiscriptiveStatistics.Max(arr) - DiscriptiveStatistics.Min(arr)) / CountInterval;
            for (int i = 0; i < m_t.Length; i++)
            {
                double u = (xn[i] - X_aver) / disp;
                double f_u = 1 / (Math.Sqrt(2 * Math.PI) * Math.Pow(Math.E, (u * u) / 2));
                m_t[i] = sum * h / disp * f_u;
            }
            return m_t;
        }

        private double[] point(double[] arr, int CountInterval)
        {
            double max = DiscriptiveStatistics.Max(arr);
            double min = DiscriptiveStatistics.Min(arr);
            double[] points = new double[CountInterval + 1];
            for (int i = 0; i < points.Length; i++)
            {
                points[i] = max - i * (max - min) / CountInterval;
            }

            return points.OrderBy(x => x).ToArray();
        }
        private double[] Get_m_e(double[] arr, double[] points, int CountInterval)
        {
            double[] m_e = new double[CountInterval];
            for (int i = 0; i < points.Length; i++)
            {
                foreach (double x in arr)
                {
                    if (points[i] < x && x < points[i + 1]) m_e[i]++;
                }
            }
            return m_e;

        }
        private double[] Get_xn(double[] arr, double[] points, int CountInterval)
        {
            double[] xn = new double[CountInterval];
            for (int i = 0; i < points.Length - 1; i++)
            {
                xn[i] = points[i + 1] - (points[i + 1] - points[i]) / 2;
            }
            return xn;
        }
        #endregion

        #region Корреляционный анализ
        private void Pair()
        {
            try
            {
                koeffPair = new DenseMatrix(colum, colum);
                for (int i = 0; i < colum; i++)
                    for (int j = i; j < colum; j++)
                        koeffPair[i, j] = koeffPair[j, i] = DiscriptiveStatistics.PairKoeff(columArray[i], columArray[j]);

                DenseMatrix t_Matrix = T_Matrix_Koeff(koeffPair);

                tbMatrix.Text = Output_R(koeffPair);
                tbMatrix.Text += "Коэффициент значим при t > 1.96 (з - значим, н - не значим)\n";

                string s = "\t";
                for (int i = 0; i < colum; i++) s += (i + 1).ToString() + "\t";
                s += "\n\n";
                for (int i = 0; i < colum; i++)
                {
                    s += (i + 1).ToString() + "\t";
                    for (int j = 0; j < colum; j++)
                        if (t_Matrix[i, j] > t_tabl)
                            s += "з\t";
                        else s += "н\t";
                    s += "\n\n";
                }
                s += "\n\n";
                tbMatrix.Text += s;
                #region Окружность
                var X_Y = new List<KeyValuePair<double, double>>()
                {
                new KeyValuePair<double,double>( 150,0),
                new KeyValuePair<double,double>( 250,40),

                new KeyValuePair<double,double>( 270,90 ),
                new KeyValuePair<double,double>( 270,190 ),

                new KeyValuePair<double, double>( 250,250 ),
                new KeyValuePair<double, double>( 150,295 ),

               new KeyValuePair<double, double>( 40,250 ),
               new KeyValuePair<double, double>( 10,190 ),

               new KeyValuePair<double, double>( 10,90 ),
               new KeyValuePair<double, double>( 40,40 )
                };
                foreach (KeyValuePair<double, double> x in X_Y)
                {
                    Ellipse l = new Ellipse();
                    l.Width = l.Height = 5;
                    l.Fill = Brushes.Red;
                    Canvas.SetTop(l, x.Value);
                    Canvas.SetLeft(l, x.Key);
                    cnvMain.Children.Add(l);
                }
                for (int i = 0; i < colum; i++)
                {
                    for (int j = i + 1; j < colum; j++)
                    {
                        if (i != j && koeffPair[i, j] >= 0.3 && koeffPair[i, j] < 0.5)
                        {
                            Line l = new Line();
                            l.X1 = X_Y[i].Key;
                            l.X2 = X_Y[j].Key;
                            l.Y1 = X_Y[i].Value;
                            l.Y2 = X_Y[j].Value;
                            l.Stroke = Brushes.Red;
                            l.StrokeThickness = 1;
                            cnvMain.Children.Add(l);
                        }
                        if (i != j && koeffPair[i, j] >= 0.5 && koeffPair[i, j] < 0.7)
                        {
                            Line l = new Line();
                            l.X1 = X_Y[i].Key;
                            l.X2 = X_Y[j].Key;
                            l.Y1 = X_Y[i].Value;
                            l.Y2 = X_Y[j].Value;
                            l.Stroke = Brushes.Yellow;
                            l.StrokeThickness = 1;
                            cnvMain.Children.Add(l);
                        }
                        if (i != j && koeffPair[i, j] >= 0.7 && koeffPair[i, j] < 1)
                        {
                            Line l = new Line();
                            l.X1 = X_Y[i].Key;
                            l.X2 = X_Y[j].Key;
                            l.Y1 = X_Y[i].Value;
                            l.Y2 = X_Y[j].Value;
                            l.Stroke = Brushes.Green;
                            l.StrokeThickness = 1;
                            cnvMain.Children.Add(l);
                        }
                    }
                }
                #endregion
            }
            catch (Exception ex) { MessageBox.Show(ex.Message, "Ошибка"); }
        }
        private void PrivateK()
        {
            try
            {
                tbMatrix2.Text = "";
                koeffPrivate = new DenseMatrix(colum, colum);
                if (koeffPair == null) throw new Exception("Вычислите парные коэффициенты!");
                for (int i = 0; i < colum; i++)
                {
                    for (int j = i + 1; j < colum; j++)
                        koeffPrivate[i, j] = koeffPrivate[j, i] = Get_AlgebralAddition(koeffPair, i, j) / Math.Sqrt(Get_AlgebralAddition(koeffPair, i, i) * Get_AlgebralAddition(koeffPair, j, j));
                    koeffPrivate[i, i] = koeffPair[i, i];
                }
                #region Окружность
                var X_Y = new List<KeyValuePair<double, double>>()
                {
                new KeyValuePair<double,double>( 150,0),
                new KeyValuePair<double,double>( 250,40),

                new KeyValuePair<double,double>( 270,90 ),
                new KeyValuePair<double,double>( 270,190 ),

                new KeyValuePair<double, double>( 250,250 ),
                new KeyValuePair<double, double>( 150,295 ),

               new KeyValuePair<double, double>( 40,250 ),
               new KeyValuePair<double, double>( 10,190 ),

               new KeyValuePair<double, double>( 10,90 ),
               new KeyValuePair<double, double>( 40,40 )
                };
                foreach (KeyValuePair<double, double> x in X_Y)
                {
                    Ellipse l = new Ellipse();
                    l.Width = l.Height = 5;
                    l.Fill = Brushes.Red;
                    Canvas.SetTop(l, x.Value);
                    Canvas.SetLeft(l, x.Key);
                    cnvMainPrivate.Children.Add(l);
                }
                for (int i = 0; i < colum; i++)
                {
                    for (int j = i + 1; j < colum; j++)
                    {
                        if (i != j && koeffPrivate[i, j] >= 0.3 && koeffPrivate[i, j] < 0.5)
                        {
                            Line l = new Line();
                            l.X1 = X_Y[i].Key;
                            l.X2 = X_Y[j].Key;
                            l.Y1 = X_Y[i].Value;
                            l.Y2 = X_Y[j].Value;
                            l.Stroke = Brushes.Red;
                            l.StrokeThickness = 1;
                            cnvMainPrivate.Children.Add(l);
                        }
                        if (i != j && koeffPrivate[i, j] >= 0.5 && koeffPrivate[i, j] < 0.7)
                        {
                            Line l = new Line();
                            l.X1 = X_Y[i].Key;
                            l.X2 = X_Y[j].Key;
                            l.Y1 = X_Y[i].Value;
                            l.Y2 = X_Y[j].Value;
                            l.Stroke = Brushes.Yellow;
                            l.StrokeThickness = 1;
                            cnvMainPrivate.Children.Add(l);
                        }
                        if (i != j && koeffPrivate[i, j] >= 0.7 && koeffPrivate[i, j] < 1)
                        {
                            Line l = new Line();
                            l.X1 = X_Y[i].Key;
                            l.X2 = X_Y[j].Key;
                            l.Y1 = X_Y[i].Value;
                            l.Y2 = X_Y[j].Value;
                            l.Stroke = Brushes.Green;
                            l.StrokeThickness = 1;
                            cnvMainPrivate.Children.Add(l);
                        }
                    }
                }
                #endregion

                tbMatrix2.Text = Output_R(koeffPrivate);
                for (int i = 0; i < colum; i++)
                {
                    for (int j = i + 1; j < colum; j++)
                    {
                        if (Math.Abs(koeffPair[i, j]) < Math.Abs(koeffPrivate[i, j])) tbMatrix2.Text += "\nКоэффициент r" + (i + 1).ToString() + "," + (j + 1).ToString() + " увеличился";
                        else tbMatrix2.Text += "\nКоэффициент r" + (i + 1).ToString() + "," + (j + 1).ToString() + " уменьшился";

                        if (koeffPrivate[i, j] * koeffPair[i, j] > 0) tbMatrix2.Text += " и изменил знак";
                        else tbMatrix2.Text += " и не изменил знак";
                    }
                }
            }
            catch (Exception ex) { MessageBox.Show(ex.Message, "Ошибка"); }
        }
        private void MultyplyK()
        {
            try
            {
                if (koeffPair == null) throw new Exception("Вычислите парные коэффициенты!");
                double det_R = koeffPair.Determinant();
                DenseMatrix R = new DenseMatrix(colum);
                double[] r_m = new double[colum];
                //tbMatrix3.Text = "Индексы детерминации для множественных коэффициентов:\n";
                //for (int i = 0; i < colum; i++) tbMatrix3.Text += Data.parametrs[i] + "\t";
                //tbMatrix3.Text += "\n";

                for (int i = 0; i < colum; i++)
                {
                    r_m[i] = Math.Sqrt(1 - (det_R / Get_AlgebralAddition(koeffPair, i, i)));
                    //tbMatrix3.Text += string.Format("{0:F2}\t", r_m[i] * r_m[i]);
                }
                //tbMatrix3.Text += "\n\nЗначимость коэффициентов:\n\n";
                double[] Fr_matrix = new double[colum];
                for (int i = 0; i < colum; i++)
                {
                    Fr_matrix[i] = r_m[i] * r_m[i] * 8 / (1 - r_m[i] * r_m[i]);
                    //tbMatrix3.Text += string.Format("{0:F2}\t", Fr_matrix[i]);
                }


                string d = "\n", r = "\n", p = "\n\n", kz = "\n", z = "\n";
                for (int i = 0; i < colum; i++)
                {
                    if (Fr_matrix[i] > F_tabl)
                    {
                        kz += string.Format("{0:F2}\n\n", Fr_matrix[i]);
                        z += string.Format("Значим\n\n");
                    }
                    else
                    {
                        kz += string.Format("{0:F2}\n\n", Fr_matrix[i]);
                        z += string.Format("Не значим\n\n");
                    }
                    r += string.Format("{0:F2}\n\n", r_m[i]);
                    d += string.Format("{0:F2}\n\n", r_m[i] * r_m[i]);
                    p += Data.parametrs[i] + "\n\n";
                }

                stpKorrDeterm.Children.Add(new TextBlock { Text = d, HorizontalAlignment = HorizontalAlignment.Center });
                stpKorrR.Children.Add(new TextBlock { Text = r, HorizontalAlignment = HorizontalAlignment.Center });
                stpKorrKoeffZN.Children.Add(new TextBlock { Text = kz, HorizontalAlignment = HorizontalAlignment.Center });
                stpKorrZN.Children.Add(new TextBlock { Text = z, HorizontalAlignment = HorizontalAlignment.Center });

                stpKorrParams.Children.Add(new TextBlock { Text = p });
            }
            catch (Exception ex) { MessageBox.Show(ex.Message, "Ошибка"); }
        }
        #region Методы
        public DenseMatrix D_Matrix(DenseMatrix Matrix)
        {
            DenseMatrix M = new DenseMatrix(colum);
            for (int i = 0; i < colum; i++)
                for (int j = 0; j < colum; j++)
                    M[i, j] = Matrix[i, j] * Matrix[i, j];
            return M;
        }
        public string Output_R(DenseMatrix Matrix)
        {
            string s = "\t";
            for (int i = 0; i < colum; i++) s += (i + 1).ToString() + "\t";
            s += "\n\n";
            for (int i = 0; i < colum; i++)
            {
                s += (i + 1).ToString() + "\t";
                for (int j = 0; j < colum; j++)
                    s += string.Format("{0:F2}\t", Matrix[i, j]);
                s += "\n\n";
            }
            s += "\n\n";
            return s;
        }
        public DenseMatrix T_Matrix_Koeff(DenseMatrix Matrix_Koeff)
        {
            DenseMatrix t_Matrix = new DenseMatrix(colum);
            for (int i = 0; i < colum; i++)
            {
                for (int j = 0; j < colum; j++)
                {
                    t_Matrix[i, j] = Math.Abs(Matrix_Koeff[i, j]) * Math.Sqrt((colum - 2) / (1 - Matrix_Koeff[i, j]));
                }
            }
            return t_Matrix;
        }
        public static DenseMatrix GetMinor(DenseMatrix A, int n, int i, int j)
        {
            DenseMatrix M = new DenseMatrix(n - 1, n - 1);

            for (int k = 0; k < n; k++)
            {
                for (int m = 0; m < n; m++)
                {
                    if (k < i && m < j) M[k, m] = A[k, m];
                    if (k > i && m > j) M[k - 1, m - 1] = A[k, m];
                    if (k > i && m < j) M[k - 1, m] = A[k, m];
                    if (k < i && m > j) M[k, m - 1] = A[k, m];
                }
            }

            return M;
        }
        private double Get_AlgebralAddition(DenseMatrix A, int i, int j)
        {
            return Math.Pow(-1, i + j) * GetMinor(A, colum, i, j).Determinant();
        }
        public DenseMatrix GetDense(double[][] A)
        {
            double[] buf = A[0];
            DenseMatrix res = new DenseMatrix(colum, buf.Length);
            for (int i = 0; i < colum; i++)
            {
                buf = A[i];
                for (int j = 0; j < buf.Length; j++) res[i, j] = buf[j];
            }
            return res;
        }
        #endregion

        #endregion

        #region Регрессионный анализ
        private void Regress()
        {
            try
            {
                double errAPR = 0;
                if (columArray == null) throw new Exception("Проведите корреляционный анализ!");
                DenseMatrix X_1 = new DenseMatrix(colum, columArray[0].Length);
                DenseMatrix x_forY = new DenseMatrix(colum, columArray[0].Length);
                DenseMatrix Y = new DenseMatrix(columArray[0].Length, 1);
                Y[0, 0] = 1.0;
                for (int i = 1; i < colum; i++)
                {
                    double[] buf = columArray[i];
                    for (int j = 0; j < columArray[0].Length; j++)
                    {
                        X_1[0, j] = 1.0;
                        x_forY[0, j] = 1.0;
                        Y[j, 0] = columArray[0][j];
                        X_1[i, j] = buf[j]; x_forY[i, j] = buf[j];
                    }
                }

                #region Коэффициенты
                DenseMatrix X_T = (DenseMatrix)X_1.Transpose();
                DenseMatrix MulMatr = (DenseMatrix)X_1.Multiply(X_T);
                DenseMatrix inverseM = (DenseMatrix)MulMatr.Inverse();
                DenseMatrix mul = (DenseMatrix)inverseM.Multiply(X_1);
                DenseMatrix A = (DenseMatrix)mul.Multiply(Y);
                B = new double[colum];
                for(int i = 0; i < colum; i++)
                {
                    B[i] = A[i, 0];
                }
                string b = string.Format("a0 = {0:F2}\n", A[0, 0]), ib = "", kz = "", z = "", y1 = "", y2 = "", y3 = "", iy = "";
                string strok = string.Format("y = {0:F2} +", A[0, 0]);
                for (int i = 1; i < A.RowCount; i++)
                {
                    b += string.Format("a{1} = {0:F2}\n", A[i, 0], i);
                    strok += string.Format(" {0:F2}*X{1} +", A[i, 0], i);
                }
                tblRegrURAVN.Text = "";
                for (int i = 0; i < strok.Length - 1; i++)
                    tblRegrURAVN.Text += strok[i];
                stpRegrToch.Children.Add(new TextBlock { Text = b, HorizontalAlignment = HorizontalAlignment.Center });
                #endregion

                #region Значимость
                DenseMatrix NewY = (DenseMatrix)X_T.Multiply(A);
                DenseMatrix Qr = (DenseMatrix)NewY.TransposeThisAndMultiply(NewY);
                double Qos = 0;
                for (int i = 0; i < columArray[0].Length; i++)
                {
                    Qos += Math.Pow(Y[i, 0] - NewY[i, 0], 2);
                }
                double S_2 = Qos / (columArray[0].Length - colum - 1);
                double S_ = Math.Sqrt(S_2);
                DenseMatrix S_b = inverseM;

                for (int i = 0; i < S_b.RowCount; i++)
                {
                    S_b[i, i] *= S_;
                    if (Math.Abs(A[i, 0]) / S_b[i, i] > t_tabl2) z += " значим\n";
                    else z += " не значим\n";
                    kz += string.Format("{0:F4}\n", Math.Abs(A[i, 0]) / S_b[i, i]);
                }
                for (int i = 0; i < colum; i++)
                {
                    ib += string.Format("{0:F2} <=A{2}<= {1:F2}\n", A[i, 0] - t_tabl2 * S_b[i, i], A[i, 0] + t_tabl2 * S_b[i, i], i);
                }

                double F = (Qr[0, 0] * (columArray[0].Length - 1 - 1)) / ((1 + 1) * Qos);
                if (F > F_tabl) tblRegrURAVN.Text += string.Format("\n Значимость равна {0:F2}. Уравнение регрессии значимо.", F);
                else tblRegrURAVN.Text += string.Format("\n Значимость равна {0:F2}. Уравнение регрессии не значимо.", F);



                double Prognoz = A[0, 0];
                DenseMatrix X0 = new DenseMatrix(1, colum);
                for (int i = 1; i < colum; i++)
                {
                    Prognoz += columArray[i][0] * A[i, 0];
                }
                X0[0, 0] = 1.0;
                for (int i = 1; i < colum; i++)
                {
                    X0[0, i] = columArray[i][0];
                }

                X0 = (DenseMatrix)X0.Transpose();

                DenseMatrix x_fory_t = (DenseMatrix)x_forY.Transpose();
                DenseMatrix MulMatr2 = (DenseMatrix)X_1.Multiply(x_fory_t);
                DenseMatrix inverseM2 = (DenseMatrix)MulMatr2.Inverse();

                DenseMatrix X_T2 = (DenseMatrix)X0.Transpose();
                DenseMatrix delta = (DenseMatrix)inverseM2.Multiply(X0);
                DenseMatrix m = (DenseMatrix)X_T2.Multiply(delta);
                for (int i = 0; i < columArray.Length; i++)
                {
                    y1 += string.Format("{0:F4}\n", Y[i, 0]);
                    y2 += string.Format("{0:F4}\n", NewY[i, 0]);
                    errAPR += Math.Abs(Y[i, 0] - NewY[i, 0]) / Y[i, 0];
                    y3 += string.Format("{0:F4}\n", Y[i, 0] - NewY[i, 0]);
                    iy += string.Format("{0:F2}<=Y^<={1:F2}\n", NewY[i, 0] - t_tabl * S_ * Math.Sqrt(Math.Abs(m[0, 0])), NewY[i, 0] + t_tabl * S_ * Math.Sqrt(Math.Abs(m[0, 0])));
                }
                errAPR /= columArray[0].Length;
                tb.Text += "\nОшибка аппроксимации: " + errAPR;
                stpRegrInterv.Children.Add(new TextBlock { Text = ib, HorizontalAlignment = HorizontalAlignment.Center });
                stpRegrIntervY_.Children.Add(new TextBlock { Text = iy, HorizontalAlignment = HorizontalAlignment.Center });
                stpRegrKoeffZN.Children.Add(new TextBlock { Text = kz, HorizontalAlignment = HorizontalAlignment.Center });
                stpRegrZN.Children.Add(new TextBlock { Text = z, HorizontalAlignment = HorizontalAlignment.Center });
                stpRegrY.Children.Add(new TextBlock { Text = y1, HorizontalAlignment = HorizontalAlignment.Center });
                stpRegrY_.Children.Add(new TextBlock { Text = y2, HorizontalAlignment = HorizontalAlignment.Center });
                stpRegrY_Y.Children.Add(new TextBlock { Text = y3, HorizontalAlignment = HorizontalAlignment.Center });
                incurracy = t_tabl * S_ * Math.Sqrt(m[0, 0] + 1);
                #endregion
            }
            catch (Exception ex) { MessageBox.Show(ex.Message, "Ошибка"); }
        }
        private void BtnRegress_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                string[] s = tbPrognozKoeff.Text.Trim().Split(' ');
                if (s.Length != colum-1) throw new Exception("Введите девять коэффициентов!");
                double[] x0 = new double[colum];
                x0[0] = 1;
                for(int i = 1; i < colum; i++)
                {
                    x0[i] = double.Parse(s[i - 1].Replace('.', ','));
                }
                double Prognoz = 0;
                for (int i = 0; i < colum; i++)
                {
                    Prognoz += x0[i] * B[i];
                }
                tblPrognoz.Text = string.Format(" Полученное знначение: Y^ = {0:F2}", Prognoz);
                tblPrognoz.Text += string.Format("\n Интервал предсказания: {0:F2}<=Y^<={1:F2}\n", Prognoz - incurracy, Prognoz + incurracy);

            }
            catch (Exception ex) { MessageBox.Show(ex.Message, "Ошибка", MessageBoxButton.OK, MessageBoxImage.Error); }
        }

        #endregion

        #region Данные
        private void BtnDATA_Click()
        {
            Data.GetCsv();
            columArraySource = Data.Array;
            columArray = new double[Data.parametrs.Count][];
            colum = Data.parametrs.Count;
            for (int i = 0; i < colum; i++)
            {
                columArray[i] = DiscriptiveStatistics.Rationing_MaxMin(columArraySource[i]);
            }
            StackPanel stp = new StackPanel();
            stp.Orientation = Orientation.Horizontal;
            for (int i = 0; i < Data.Array.Length; i++)
            {
                double[] buf = Data.Array[i];
                TextBlock lv = new TextBlock();
                lv.Width = 100;
                lv.Height = Data.Array[0].Length * 18;
                lv.Text += Data.parametrs[i] + "\n";
                for (int j = 0; j < buf.Length; j++)
                    lv.Text += buf[j] + "\n";
                stp.Children.Add(lv);
            }
            scv.Content = stp;

            string s = "";

            for (int j = 0; j < columArray[0].Length; j++)
            {
                for (int i = 0; i < colum; i++)
                {
                    s += columArray[i][j] + " ";
                }
                s += "\n";
            }
            tb.Text = s;
        }
        #endregion

        #region Описательная статистика
        private void BtnDiscrStat_Click(object sender, RoutedEventArgs e)
        {
            stpDiscrStat.Children.Clear();
            try
            {
                if (rbNorm.IsChecked == true)
                {
                    output(columArray);
                    tblDiscrStat.Text = "\n\nСреднее\n\nСумма\n\nСтандартная ошибка\n\nМедиана\n\nМода\n\nСтандартное отклонение\n\nДисперсия\n\nЭксцесс\n\nАсимметричность\n\nИнтервал\n\nМинимум\n\nМаксимум\n\nСчет";
                }
                else if (rbNotNorm.IsChecked == true)
                {
                    output(columArraySource);
                    tblDiscrStat.Text = "\n\nСреднее\n\nСумма\n\nСтандартная ошибка\n\nМедиана\n\nМода\n\nСтандартное отклонение\n\nДисперсия\n\nЭксцесс\n\nАсимметричность\n\nИнтервал\n\nМинимум\n\nМаксимум\n\nСчет";
                }
                else throw new Exception("Выберите вид данных!");
            }
            catch (Exception ex) { MessageBox.Show(ex.Message, "Ошибка"); }
        }

        private void output(double[][] arr)
        {
            for (int i = 0; i < colum; i++)
            {
                string s = "";
                string t = "";
                s += (i + 1).ToString() + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Average(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Amount(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.StandartError(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Median(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Fashion(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.StandartDeviation(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Dispersion(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Excess(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Asymmetry(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Interval(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Min(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", DiscriptiveStatistics.Max(arr[i])) + "\n\n";
                t += "\t\n\n\n";
                s += string.Format("{0:F2}", arr[i].Length) + "\n\n";
                t += "\t\n\n\n";
                stpDiscrStat.Children.Add(new TextBlock { Text = s });
                stpDiscrStat.Children.Add(new TextBlock { Text = t });
            }
        }


        #endregion

        private void err()
        {
            double srErr = 0;
            
            tb.Text = "";
            for (int i = 0; i < colum; i++)
            {
                double t = 2.5;
                double[] buf = columArray[i];
                double disp = DiscriptiveStatistics.Dispersion(buf);
                srErr += DiscriptiveStatistics.StandartError(buf);
                double err = Math.Sqrt(disp / buf.Length);
                double err2 = t * err;
                double n = t * t * disp * buf.Length / (err2 * err2 * buf.Length + t * t * disp);
                tb.Text += i + "\tn = " + n + "\tсред.ош. = " + err + "\tпредел.ош. = " + err2 + "\n";
            }
            srErr /= colum;
            tb.Text += "\nСредняя ошибка: " + srErr;

        }
    }
}
