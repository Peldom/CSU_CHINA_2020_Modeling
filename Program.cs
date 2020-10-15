using System;
using System.IO;
using System.Text;

namespace Alga_Dynamics
{
    class Program_Main
    {
        static void Main(string[] args)
        {
            alginate_ball b1 = new alginate_ball(0.47156, 1,0.5,10000);
            b1.run_simulate(80);
        }
    }
    class alginate_ball
    {
        public int ball_num { get; set; }
        public double env_cd { get; set; }// mg/cm^3
        public double env_v { get; set; }

        public double ball_radius { get; set; } //cm
        public static double aq_radius = 1.0 * Math.Pow(10, -4); //cm
        public double aq_density = Math.Pow(10, 7) / 6.0;//cell/cm^3
        public int shell_num{ get; set; }
        public double shell_radius_unit { get; set; }
        public double[] Cd_shells;
        public aquetilis_cell[] aq_shell;
        public double[] V_shells;
        public alginate_ball(double radius,int ball_num,double env_cd,double env_v)
        {
            ball_radius = radius;
            this.ball_num = ball_num;
            this.env_cd = env_cd;
            this.env_v = env_v;
            shell_num = (int)(ball_radius / aq_radius)/100;
            shell_radius_unit = ball_radius / shell_num;
            double[] shell_aq_num = new double[shell_num];
            V_shells = new double[shell_num];
            Cd_shells = new double[shell_num];
            aq_shell = new aquetilis_cell[shell_num];//null array yet
            for (int i = 0; i < shell_num; i++)
            {
                V_shells[i] = Add_Math.v_ball(shell_radius_unit * (i + 1)) - Add_Math.v_ball(shell_radius_unit * i);
                shell_aq_num[i] = V_shells[i] * aq_density;
                Cd_shells[i] = 0;//default                
                aq_shell[i] = new aquetilis_cell(shell_aq_num[i]);
            }
            //foreach (var a in shell_aq_num)
            //{
            //    Console.WriteLine(a);
            //}
            //Console.WriteLine("Shell Num = " + shell_num);
            //Console.WriteLine("aquetilis cell number per ball = " + Add_Math.v_ball(radius) * aq_density);
        }     
        public double Shell_Cd_infiltrate(double Cd_p1,int shell_index,double Cd_p2)
        {
            //assuming Cd_p2 V > Cd_p1_V
            //attention! it could also be negative
            //fisk law
            double contact_surface =Add_Math.s_circle((1 + shell_index) * shell_radius_unit);
            return contact_surface*(Cd_p2-Cd_p1)*0.0015;
            //lack of 0 protect
        } 
        public double sum_env_infiltrate_ball(double env_cd,int ball_num,double ball_surface_cd)
        {
            //undifined
            //attention! it could also be negative
            //fisk law
            double contact_surface = Add_Math.s_circle(ball_radius) * ball_num;//sum surface
            double Cd_trans = (env_cd - ball_surface_cd) * 0.0025;
            return Cd_trans * contact_surface;
            //return 0.01;
        }
        public void run_simulate(int ts)  
        {
            for (int i = 0; i < ts; i++)
            {
                //ball absorb Cd from environment
                double sum_env_cd_num_to_ball = sum_env_infiltrate_ball(env_cd, ball_num, Cd_shells[shell_num - 1]);//mg
                env_cd -= sum_env_cd_num_to_ball/env_v;//mg/cm3
                Cd_shells[shell_num - 1] += sum_env_cd_num_to_ball /ball_num/V_shells[shell_num-1];//outest,mg/cm3              

                //only for writecsv
                double[] aq_num_shell = new double[shell_num];
                double[] aq_Cd_in = new double[shell_num];
                //double[] aq_death = new double[shell_num];
                double sum_shell_cd=0;
                double sum_cell_cd=0;

                for (int shell_index = shell_num-1; 0<shell_index; shell_index--)
                {
                    //Cd infiltrate between shells
                    //defaut positive
                    //from outside to inside;might be a huge logic bug on science! (Direction: only one)
                    double Cd_auto_infl = Shell_Cd_infiltrate(Cd_shells[shell_index-1], shell_index, Cd_shells[shell_index]);                  
                    Cd_shells[shell_index] -= Cd_auto_infl/V_shells[shell_index];
                    Cd_shells[shell_index - 1] += Cd_auto_infl/ V_shells[shell_index-1];

                    //aquetilis absorb Cd
                //    double Cd_aq_absorb = aq_shell[shell_index].transfer(Cd_shells[shell_index]);
                   // if(Cd_shells[shell_index] -Cd_aq_absorb * aq_shell[shell_index].num > 0)
                    //{
                  //      Cd_shells[shell_index] -= Cd_aq_absorb * aq_shell[shell_index].num;
                    //    aq_shell[shell_index].Cd_in += Cd_aq_absorb;
                    //}
                 

                    //aquetilis number change
                  //  aq_shell[shell_index].re_gen(Cd_shells[shell_index], aq_shell[shell_index].num);

                    //aquetilis Cd death feedback
                    //defaut posi tive
                    //double Cd_aq_death_release = aq_shell[shell_index].death_release();
                    //Cd_shells[shell_index] += Cd_aq_death_release;

                    //csv array write
                    aq_num_shell[shell_index] = aq_shell[shell_index].num;
                    aq_Cd_in[shell_index] = aq_shell[shell_index].Cd_in;// aq_shell[shell_index].num;
                    //aq_death[shell_index] = aq_shell[shell_index].death_rate;
                    sum_cell_cd += aq_shell[shell_index].Cd_in* aq_shell[shell_index].num* Add_Math.v_ball(aq_radius);//mg
                    sum_shell_cd += Cd_shells[shell_index]*V_shells[shell_index];
                }
                //aq_num_shell[0] = aq_shell[0].num;
                aq_num_shell[0] = aq_shell[0].num;
                aq_Cd_in[0] = aq_shell[0].Cd_in;// aq_shell[0].num;
                aq_shell[0].re_gen(Cd_shells[0], aq_shell[0].num);
                sum_cell_cd += aq_shell[0].Cd_in*aq_shell[0].num*Add_Math.v_ball(aq_radius);
                sum_shell_cd += Cd_shells[0]*V_shells[0];
                Console.WriteLine("{0} Env_Cd:{1:f4}mg\tsum_cell_cd:{2:f4}mg\tsum_shell_cd:{3:f4}mg\tshell_outest_cd:{4:f2}mg/cm3\toutest_aq_cd:{5:f2}mg/cm3", i, env_cd * env_v, sum_cell_cd, sum_shell_cd, Cd_shells[shell_num - 1], aq_shell[shell_num - 1].Cd_in);// * aq_shell[shell_num - 1].num);
                filewrite.WriteCSV(Cd_shells, shell_num);
            }
        }
    }
    class aquetilis_cell
    {    
        public double Cd_in { get; set; }//of each cell,not cell group
        public double death_rate { get; set; }
 //       public double growth_rate { get; set; }
        public double num { get; set; }
        public aquetilis_cell(double previous_num)
        {
            Cd_in =0;
            num = previous_num;
        }
        public void re_gen(double Cd_out, double previous_num)
        {
            Cd_in += transfer(Cd_out) / previous_num/Add_Math.v_ball(alginate_ball.aq_radius);
            death_rate = 0.05 + death_change();   
 //           growth_rate = 0.1 + growth_change();
            num = previous_num * (1-death_rate);
        }
        private double death_change()
        {
            //todo
            return Cd_in*0.01;
        }
        //private double growth_change()
        //{
        //    //todo
        //    return -Cd_in*0.02;
        //}
        public double transfer(double Cd_out)//of each cell
        {//assert transfer>0

            if (Cd_in < Cd_out * 100) return 0.0000000001;
            else return 0;    
        }
        public double death_release()
        {
            return death_rate*num * 0.2 * Cd_in;
            //return 0;
        }
    }
    class Add_Math
    {
        public static double v_ball(double radius)
        {
            return 4/3*Math.PI * Math.Pow(radius, 3.0);//Π*r^3*4/3
        }
        public static double s_circle(double radius)
        {
            return Math.PI * Math.Pow(radius, 2.0) * 4;//Π*r^2*4
        }
    }
    class filewrite
    {
        public static void WriteCSV(double[] shell_aq_num, int shell_num)
        {
            string path = @"D:\Remain\iGEM\model\种群\CSVinfoalltest.csv";
            if (!File.Exists(path)) File.Create(path).Close();
            StreamWriter sw = new StreamWriter(path, true, Encoding.UTF8);
            for (int i = 0; i < shell_num; i++)
            {
                sw.Write(shell_aq_num[i] + ",");
            }
            sw.Write("\r\n");
            //sw.Write(sch_len + ",");
            //sw.Write(data_ds[1] + "\r\n");
            sw.Flush();
            sw.Close();
        }
    }
} 