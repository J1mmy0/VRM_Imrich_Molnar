using System;
using System.Linq;
using UnityEngine;
using MathNet.Numerics.LinearAlgebra;

public class Core : MonoBehaviour
{
    public static class GlobalVariables
    {
        public static float[] Q_Target = new float[2];
        public static float[] Q_Actual = new float[2];
        public static bool[] Q_In_Pos = new bool[2];
    }

    private double EE_x, EE_y;
    private double Q_0, Q_1;
    public GameObject Obj_Joint_0, Obj_Joint_1;
    public GameObject Target_EE;
    public bool In_Position;
    public bool Execute;

    private static readonly double[] Link_length = new double[2] {1.0, 0.5};
    private static readonly double[,] Joint_limit = new double[2, 2] {{-180.0, 180.0 }, 
                                                                      {-180.0, 180.0 }};

    private int state_id = 0;

    void Start()
    {
        Obj_Joint_0.transform.localEulerAngles = new Vector3(0.0f, 0.0f, 0.0f);
        Obj_Joint_1.transform.localEulerAngles = new Vector3(0.0f, 0.0f, 0.0f);

        (EE_x, EE_y) = FK((-1.0) * Degree_To_Radian(Obj_Joint_0.transform.localEulerAngles.z), (-1.0) * Degree_To_Radian(Obj_Joint_1.transform.localEulerAngles.z));
        Target_EE.transform.localPosition = new Vector3((float)EE_x, 0.3f, (float)EE_y);

        In_Position = true; Execute = false; 
    }

    void Update()
    {
        switch (state_id)
        {
            case 0:
                {
                    In_Position = true;

                    if (Execute == true)
                    {
                        EE_x = Target_EE.transform.localPosition.x;
                        EE_y = Target_EE.transform.localPosition.z;

                        state_id = 1;
                    }
                }
                break;

            case 1:
                {
                    Execute = false; In_Position = false;
                    GlobalVariables.Q_In_Pos[0] = false; GlobalVariables.Q_In_Pos[1] = false;

                    (Q_0, Q_1) = IK((double)EE_x, (double)EE_y, (-1.0) * Degree_To_Radian(Obj_Joint_0.transform.localEulerAngles.z), (-1.0) * Degree_To_Radian(Obj_Joint_1.transform.localEulerAngles.z), 
                                    0.0001f, 10000);

                    // Ensure joint angles stay within the limit of -180 to 180 degrees
                    Q_0 = ClampAngle(Q_0, -180, 180);
                    Q_1 = ClampAngle(Q_1, -180, 180);

                    Obj_Joint_0.transform.localEulerAngles = new Vector3(0.0f, 0.0f, (float)((-1) * Radian_To_Degree(Q_0)));
                    Obj_Joint_1.transform.localEulerAngles = new Vector3(0.0f, 0.0f, (float)((-1) * Radian_To_Degree(Q_1)));

                    state_id = 0;
                }
                break;

            case 2:
                {
                    if (GlobalVariables.Q_In_Pos[0] == true && GlobalVariables.Q_In_Pos[1] == true)
                    {
                        state_id = 0;
                    }
                }
                break;
        }
    }

    public static (double, double) FK(double Q_0, double Q_1)
    {
        return (Link_length[0] * Math.Cos(Q_0) + Link_length[1] * Math.Cos(Q_0 + Q_1),
                Link_length[0] * Math.Sin(Q_0) + Link_length[1] * Math.Sin(Q_0 + Q_1));
    }

    public static (double, double) IK(double TCP_x, double TCP_y, double Q_0_actual, double Q_1_actual, double accuracy, int N)
    {
        double[] TCP_target = { TCP_x, TCP_y };
        double[] Q_actual = { Degree_To_Radian(Q_0_actual), Degree_To_Radian(Q_1_actual) };

        (double x, double y) = FK(Q_actual[0], Q_actual[1]);

        double[] TCP_actual = { x, y };
        double[] e_p = {TCP_target[0] - TCP_actual[0],
                        TCP_target[1] - TCP_actual[1]};

        foreach (var _ in Enumerable.Range(0, N))
        {
            Matrix<double> J = Get_Jacobian(Q_actual);
            Vector<double> v_e = Vector<double>.Build.Dense(new double[] { e_p[0], e_p[1] });
            Vector<double> Q_i = J.PseudoInverse() * v_e;

            Q_actual[0] += Q_i[0];
            Q_actual[1] += Q_i[1];

            (TCP_actual[0], TCP_actual[1]) = FK(Q_actual[0], Q_actual[1]);

            e_p[0] = TCP_target[0] - TCP_actual[0];
            e_p[1] = TCP_target[1] - TCP_actual[1];

            if (Euclidean_Norm(e_p) < accuracy)
            {
                break;
            }
        }

        return (Q_actual[0], Q_actual[1]);
    }

    public static Matrix<double> Get_Jacobian(double[] Q)
    {
        Matrix<double> J = Matrix<double>.Build.DenseOfArray(new double[,]
        {{ 0.0, 0.0}, {0.0, 0.0} });

        J[0, 0] = (-1) * Link_length[0] * Math.Sin(Q[0]) + ((-1) * Link_length[1] * Math.Sin(Q[0] + Q[1]));
        J[0, 1] = (-1) * Link_length[1] * Math.Sin(Q[0] + Q[1]);
        J[1, 0] = Link_length[0] * Math.Cos(Q[0]) + (Link_length[1] * Math.Cos(Q[0] + Q[1]));
        J[1, 1] = Link_length[1] * Math.Cos(Q[0] + Q[1]);

        return J;
    }

    public static double Degree_To_Radian(double x)
    {
        return x * (Math.PI / 180.0);
    }

    public static double Radian_To_Degree(double x)
    {
        return x * (180.0 / Math.PI);
    }

    public static double Euclidean_Norm(double[] x)
    {
        double x_sum = 0;
        foreach (double x_i in x)
        {
            x_sum += Math.Pow(x_i, 2);
        }

        return Math.Sqrt(x_sum);
    }

    // Function to clamp angle between min and max
    public static double ClampAngle(double angle, double min, double max)
    {
        if (angle < min)
            angle += 360;
        else if (angle > max)
            angle -= 360;

        return Mathf.Clamp((float)angle, (float)min, (float)max);
    }
}