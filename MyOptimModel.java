import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.experimental.ParallelComputer;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.UnknownObjectException;

public class MyOptimModel {
	// Constant
	static final int omega = 2; // 1 => 1-min, 2=>max,  else => 1-avg
	static final int mode = 0; // 1 => A and B, 0 => othw
	static final double cx = 7000.0;
	static final double cy = 7000.0;
	static final int nofgat =1;
	static final int nofcli = 50;
	static final Double Sigma = 0.5; // how much gateways % is not possible to exceed
	static final double threshold = 0.50; // threshold for Pr(i)
	static final double beta = 0.66; // threshold for H : 66% for large range, 95% for short range
	static final int BW = 125; // in kHz
	static final int payload = 59; // in byte
	static final double lambda = 0.0071; // # packet / second
	static final double Ptx = 14.0; // in dBm
	static final Integer[] SF = { 7, 8, 9, 10, 11, 12 };
	static final Double[] SNR = { -6.0, -9.0, -12.0, -15.0, -17.5, -20.0 }; // in dB
	static final Double[][] SINR = { { 6.0, -16.0, -18.0, -19.0, -19.0, -20.0 },
			{ -24.0, 6.0, -20.0, -22.0, -22.0, -22.0 }, { -27.0, -27.0, 6.0, -23.0, -25.0, -25.0 },
			{ -30.0, -30.0, -30.0, 6.0, -26.0, -28.0 }, { -33.0, -33.0, -33.0, -33.0, 6.0, -29.0 },
			{ -36.0, -36.0, -36.0, -36.0, -36.0, 6.0 } }; // capture effect power thresholds per SF[SF] for
															// non-orthogonal transmissions
	static final Double[] ToA = { 102.7, 184.8, 328.7, 616.5, 1315.0, 2466.0 }; // in ms

	// SETS
	private Set<Node> sensors;
	private Set<Node> gateways;

	static Double Havg;

	// Constants
	static final double eps = 1.0e-7;
	static final double M = 1.0e5; // big number

	// ILP
	private IloCplex cplex;

	// Variables

	private HashMap<Node, IloNumVar> k; // Map between name gateway and if the gateway is deployed or not [1/0]
	private Map<Node, Map<Integer, IloNumVar>> y; // Map between <node : <SF : node allocated with SF f>>
	private Map<Node, Map<Node, Map<Integer, IloNumVar>>> j; // Variable Jxkyif
	private HashMap<Node, Map<Node, Map<Integer, IloNumVar>>> a;// Variable A ij f
	private HashMap<Node, Map<Node, Map<Integer, Map<Integer, IloNumVar>>>> b; // Variable B ij ff'
	private HashMap<Node, Map<Node, Map<Node, Map<Integer, Map<Integer, IloNumVar>>>>> w; // Varibale WJ
	private HashMap<Node, Map<Integer, IloNumVar>> s; // Variable Sif
	private HashMap<Node, Map<Integer, Map<Node, Map<Integer, IloNumVar>>>> r; // Variable RyfSif

	// Parameters
	private HashMap<String, HashMap<Integer, Double>> Hs; // H for each Gateway/Sensor and f pair.
	private HashMap<String, Double> Prx_i_k; // PRX for the sensor i allocated to node k
	private HashMap<String, HashMap<String, HashMap<Integer, Integer>>> Cf; // Collision same SF for node i allocated
	private HashMap<String, HashMap<String, HashMap<String, Integer>>> If; // to k1 and j allocated to k2 at sf F.
	private HashMap<String, HashMap<String, HashMap<Integer, Double>>> Hss;
	private HashMap<Node, Map<Node, Map<Integer, IloNumVar>>> z;

	public MyOptimModel() throws IloException {
		// Sets
		sensors = new HashSet<Node>();
		gateways = new HashSet<Node>();
		cplex = new IloCplex();

		// Variables
		y = new HashMap<Node, Map<Integer, IloNumVar>>(); // Classic for the opt function
		k = new HashMap<Node, IloNumVar>(); // For understanding if k gateway is activated or not
		j = new HashMap<Node, Map<Node, Map<Integer, IloNumVar>>>(); // Variable JxkYif
		w = new HashMap<Node, Map<Node, Map<Node, Map<Integer, Map<Integer, IloNumVar>>>>>(); // Variable WJxkYifYjf
		a = new HashMap<Node, Map<Node, Map<Integer, IloNumVar>>>(); // Variable A ij f
		b = new HashMap<Node, Map<Node, Map<Integer, Map<Integer, IloNumVar>>>>(); // Variable B ij ff'
		s = new HashMap<Node, Map<Integer, IloNumVar>>(); // GET MAX H
		z = new HashMap<Node, Map<Node, Map<Integer, IloNumVar>>>();
		r = new HashMap<Node, Map<Integer, Map<Node, Map<Integer, IloNumVar>>>>();

		// Parameters
		Hs = new HashMap<String, HashMap<Integer, Double>>(); // H for each n-uple Sensor-Gateway-SNR
		Hss = new HashMap<String, HashMap<String, HashMap<Integer, Double>>>();
		// Prx_i_k = new HashMap<String,Double>(); //For Prx for the node i to the
		// gateway k
		Cf = new HashMap<String, HashMap<String, HashMap<Integer, Integer>>>();
		If = new HashMap<String, HashMap<String, HashMap<String, Integer>>>();

	}

	public void setSensors(double x, double y, int nbs) { // Generate randomly sensors in x y.
		sensors = new HashSet<Node>();
		// Generate sensor's position randomly on the area within radius xmax/2 and H >
		// beta
		Node g = new Node(x / 2, y / 2, x / 2, "Gtemp"); // create a node temp and assign it x/2 y/2 and x/2 as ptx
		for (int i = 0; i < nbs; i++) {
			double xtemp = Math.random() * x;
			double ytemp = Math.random() * y;
			double x2 = x / 2; // era y/2
			double y2 = y / 2; // era x/2
			double dist = Math.sqrt((xtemp - x2) * (xtemp - x2) + (ytemp - y2) * (ytemp - y2));
			Node temp = new Node(xtemp, ytemp, Ptx, "S" + i);
			while (dist > (x / 2) || temp.getH(g, SNR[5]) < beta) { // snr 5 is the highest
				// System.out.println("i = "+i+" "+dist+ " H(SF12) = "+temp.getSuccessProba(g,
				// SNR[5]));
				xtemp = Math.random() * x;
				ytemp = Math.random() * y;
				dist = Math.sqrt((xtemp - x2) * (xtemp - x2) + (ytemp - y2) * (ytemp - y2));
				temp = null;
				temp = new Node(xtemp, ytemp, Ptx, "S" + i);
			}
			sensors.add(new Node(xtemp, ytemp, Ptx, "S" + i));
		}
	}

	public void setGatewaysR(double x, double y, int nbg) { // Gateways are generated randomly in the area x y
		gateways = new HashSet<Node>();
		int id = 0;

		for (int i = 0; i < nbg; i++) {

			double xtemp = 0 + (double) Math.random() * (x - 0); // Generate a number for x between 0 and 2x
			double ytemp = 0 + (double) Math.random() * (y - 0); // same for y
			gateways.add(new Node(xtemp, ytemp, x / 2, "G" + id));
			id++;
		}
	}

	public void ComputeH() {

		for (Node sensor : sensors) {
			HashMap<Integer, Double> prov = new HashMap<Integer, Double>();
			for (Node gateway : gateways) {
				for (int sf = 7; sf <= 12; sf++) {
					prov.put(sf, sensor.getH(gateway, SNR[sf - 7]));
					System.out.println(
							sensor.getName() + gateway.getName() + sf + " " + sensor.getH(gateway, SNR[sf - 7]));
				}
				Hs.put(sensor.getName() + gateway.getName(), prov);
			}

		}
	}

	public void ComputeHss() {
		//Havg = (double) 0;
		//int count = 0;
		for (Node sensor : sensors) {
			HashMap<String, HashMap<Integer, Double>> prov1 = new HashMap<String, HashMap<Integer, Double>>();
			for (Node gateway : gateways) {
				HashMap<Integer, Double> prov = new HashMap<Integer, Double>();
				for (int sf = 7; sf <= 12; sf++) {
					prov.put(sf, sensor.getH(gateway, SNR[sf - 7]));
					
					//Havg += sensor.getH(gateway, SNR[sf - 7]);
					//count++;
					 System.out.println(sensor.getName() + gateway.getName() + sf +" "+
					 sensor.getH(gateway, SNR[sf - 7]));
				}
				prov1.put(gateway.getName(), prov);
			}
			Hss.put(sensor.getName(), prov1);
		}
		//Havg = Havg / count;
		// System.out.println("Haveg = "+Havg+"\n");
	}

	public void PrintH_totxt() throws IOException {
		String f = "./H.txt/";
		BufferedWriter wout = new BufferedWriter(new FileWriter(f));
		wout.write("Calcolo H\n");
		for (String key1 : Hs.keySet()) {
			wout.write("\n Sensore: " + key1 + ":\n");
			for (int key2 : Hs.get(key1).keySet()) {
				wout.write("SF : " + key2 + "= " + Hs.get(key1).get(key2) + " <<>> ");
			}
		}
		wout.close();
	}

	public void printHss_totxt() throws IOException {
		String f = "./Hss.txt/";
		BufferedWriter wout = new BufferedWriter(new FileWriter(f));
		wout.write("Calcolo Hss\n");
		for (String key1 : Hss.keySet()) {
			for (String key2 : Hss.get(key1).keySet()) {
				for (int key3 : Hss.get(key1).get(key2).keySet()) {
					wout.write("\n Sensor" + key1 + " and Gateway : " + key2 + "; Sf: " + key3 + " "
							+ Hss.get(key1).get(key2).get(key3) + "\n");
				}
			}

		}
		wout.close();
	}

	public void PrintPrx_totxt() throws IOException {
		String f = "./Prx.txt/";
		BufferedWriter wout = new BufferedWriter(new FileWriter(f));
		wout.write("Calcolo Prx\n");
		for (String key1 : Prx_i_k.keySet()) {
			wout.write("\n Sensor and Gateway : " + key1 + "; Value: " + Prx_i_k.get(key1) + "\n");
		}
		wout.close();
	}

	public void ComputeCf() {
		// Cf = new HashMap<String, HashMap<String, HashMap <Integer, Integer>>>
		for (Node s1 : sensors) {

			for (Node s2 : sensors) {
				HashMap<String, HashMap<Integer, Integer>> cfpg1g2 = new HashMap<String, HashMap<Integer, Integer>>();
				for (Node k : gateways) {
					HashMap<Integer, Integer> cfp = new HashMap<Integer, Integer>();
					for (Integer sf : SF) {

						if (!(s1.equals(s2))) // Sensor must be different but gateway can be the same
						{
							if (Node.invDb(s1.getReceivedPower(k)) / Node.invDb(s2.getReceivedPower(k)) - 4 <= eps) {
								cfp.put(sf, 1);
							} else {
								cfp.put(sf, 0);
							}
						}

						cfpg1g2.put(k.getName(), cfp);

						// Cf.put(s1.getName()+s2.getName(), cfpg1g2);

					}

				}

				Cf.put(s1.getName() + s2.getName(), cfpg1g2);

			}

		}
	}

	public void ComputeIf() {
		for (Node s1 : sensors) {
			for (Node s2 : sensors) {
				HashMap<String, HashMap<String, Integer>> cfpg1g2 = new HashMap<String, HashMap<String, Integer>>();
				for (Node k : gateways) {
					HashMap<String, Integer> cfp = new HashMap<String, Integer>();
					for (Integer sf1 : SF) {
						for (Integer sf2 : SF) {
							if (!(s1.equals(s2))) // Sensor must be different but gateway can be the same
							{
								if (!sf1.equals(sf2)) {

									if (s1.getReceivedPower(k) - s2.getReceivedPower(k) <= SINR[sf1 - 7][sf2 - 7]) {
										if (true) {
											cfp.put(sf1.toString() + "_" + sf2.toString(), 1);
										}

									} else {
										if (true)
											cfp.put(sf1.toString() + "_" + sf2.toString(), 0);
									}
								}
							}
						}

					}

					cfpg1g2.put(k.getName(), cfp);

				}

				If.put(s1.getName() + s2.getName(), cfpg1g2);

			}

		}
	}

	public void PrintIF_totxt() throws IOException {
		String f = "./IF.txt/";
		BufferedWriter wout = new BufferedWriter(new FileWriter(f));
		wout.write("Computation of IF \n");
		for (String k1 : If.keySet()) {
			for (String k2 : If.get(k1).keySet()) {
				for (String k3 : If.get(k1).get(k2).keySet()) {
					wout.write("\n");
					wout.write(k1 + "_" + k2 + ":" + "F: =  " + k3 + " " + "=" + If.get(k1).get(k2).get(k3));
				}
			}
		}
		wout.close();
	}

	public void PrintCF_totxt() throws IOException {
		String f = "./CF.txt/";
		BufferedWriter wout = new BufferedWriter(new FileWriter(f));
		wout.write("Computation of CF \n");
		for (String k1 : Cf.keySet()) {
			for (String k2 : Cf.get(k1).keySet()) {
				for (int k3 : Cf.get(k1).get(k2).keySet()) {
					wout.write("\n");
					wout.write(k1 + "_" + k2 + ":" + "F: " + k3 + " " + Cf.get(k1).get(k2).get(k3));
				}
			}
		}
		wout.close();
	}

	public void coh() {
		for (Node s1 : sensors) {
			for (Node s2 : sensors) {
				if (!s1.equals(s2)) {
					for (Node k : gateways) {

						double cik = Node.invDb(s1.getReceivedPower(k)) / Node.invDb(s2.getReceivedPower(k));
						String s;
						if (cik <= 6)
							s = "Interfer";
						else
							s = "Not intefer";
						System.out.print("Power betw Sensor " + s1.getName() + " and Sensor " + s2.getName() + " is "
								+ cik + " so they will " + s + " at gateway " + k.getName() + "\n");

					}
				}
			}
		}

	}

	public static void main(String[] args) throws IloException, IOException {

		MyOptimModel model = new MyOptimModel();
		model.setGatewaysR(cx, cy, nofgat);
		model.setSensors(cx, cy, nofcli);
		// model.ComputeH();
		model.ComputeHss();
		model.printHss_totxt();
		// model.PrintH_totxt();
		// model.computePrx();
		// model.PrintPrx_totxt();
		model.ComputeCf();
		model.PrintCF_totxt();
		model.ComputeIf();
		model.PrintIF_totxt();
		model.setModel();
		// model.coh();
		// model.prova();

		// Print sensors
		/*
		 * graphic example = new graphic("Gateways and Clients", model.sensors,
		 * model.gateways); example.setSize(800, 400);
		 * example.setLocationRelativeTo(null); example.setVisible(true);
		 */
	}
	// MISS PRINT STAT TOPO

	/*
	 * 
	 * SET COllisions
	 * 
	 * 
	 */
	public boolean IntraSF(Node i, Node j, int sf) {
		boolean temp = !i.equals(j) && y.get(i).containsKey(sf) && y.get(j).containsKey(sf);
		boolean gtws = true;
		for (Node g : gateways)
			gtws = gtws && Node.invDb(i.getReceivedPower(g)) / Node.invDb(j.getReceivedPower(g)) - 4 <= eps;
		return temp && gtws;
	}
	/*
	 * 
	 * 
	 * LINEAR MODEL
	 */

	public double getH_obj(Node n, int i) {

		String n1 = n.getName();
		if (i == 1) // min
		{
			double Thold = 10;

			for (String s1 : Hss.get(n1).keySet()) {
				for (int s2 : Hss.get(n1).get(s1).keySet()) {
					if (Hss.get(n1).get(s1).get(s2) < Thold)
						Thold = Hss.get(n1).get(s1).get(s2);
				}
			}
			//System.out.println(Thold);
			return (Thold); // mod here hh
		} else if(i ==2)
		{
			double Thold = 0;
			for (String s1 : Hss.get(n1).keySet()) {
				for (int s2 : Hss.get(n1).get(s1).keySet()) {
					if (Hss.get(n1).get(s1).get(s2) > Thold)
						Thold = Hss.get(n1).get(s1).get(s2);
				}
			}
			//System.out.println(Thold);
			return (1-Thold);
		}
		
		else // avg
		{
			double add = 0;
			int count = 0;
			for (String s1 : Hss.get(n1).keySet()) {
				for (int s2 : Hss.get(n1).get(s1).keySet()) {
					add = add + Hss.get(n1).get(s1).get(s2);
					count++;
				}

			}
			return (Double) ( (add / count));
		}

	}

	public void setObjective() throws IloException {
		IloLinearNumExpr obj = cplex.linearNumExpr();
		// lambda <= (1 - \sum_f y^f_i) + \sum_f H y^f_i

		for (Node s : sensors) { // Setto il medio
			for (int f : y.get(s).keySet()) {
				obj.addTerm(getH_obj(s, omega), y.get(s).get(f));
				//obj.addTerm(1, y.get(s).get(f));
				// PROVARE SENZA W
				// obj.addTerm(1, r.get(s).get(f).get(s).get(f));//Modifica per allocare sf
				// basse vicino ai gateway
				// System.out.println(getH_obj(s,2) + "\n");
			}
		}
		cplex.addMaximize(obj);
	}

	void setModel() throws IloException, IOException {
		cplex.setParam(IloCplex.Param.TimeLimit, 3600);
		cplex.setOut(null);
		cplex.setWarning(null);
		setLP();
		setObjective();
		setCns(mode);
		solve("solved.txt");

	}

	public void setCns(int mode) throws IloException {
		// Only one Spfactor foreach sensor
		for (Node i : sensors) {
			// sum_sf y^sf_i <= 1
			IloLinearNumExpr sfCstr = cplex.linearNumExpr();
			for (int sf : y.get(i).keySet()) {
				sfCstr.addTerm(1.0, y.get(i).get(sf));
			}
			cplex.addLe(sfCstr, 1);
		}

		// Gateway Constraint
		IloLinearNumExpr gatcns = cplex.linearNumExpr();
		for (Node gat : gateways) {
			gatcns.addTerm(1, k.get(gat));
			// System.out.println("Inserisco " + gat.getName()+"\n");
		}
		cplex.addLe(gatcns, Math.ceil(k.size() * Sigma));

		// Linearizzazione di Xk Yi constraint

		for (Node i : sensors) {
			for (Node gat : gateways) {
				for (int sf : SF) {
					cplex.addLe(j.get(i).get(gat).get(sf), y.get(i).get(sf)); // jx_ky_i^f <= y_i^f
					cplex.addLe(j.get(i).get(gat).get(sf), k.get(gat));// jx_ky_i^f <= X_k
					IloLinearNumExpr lincns = cplex.linearNumExpr();
					IloNumVar n = cplex.numVar(1, 1);
					lincns.addTerm(1, k.get(gat));
					lincns.addTerm(1, y.get(i).get(sf));
					lincns.addTerm(-1, n);
					cplex.addGe(j.get(i).get(gat).get(sf), lincns);
				}
			}
		}

		// Zikf leq constraint
		for (Node i : sensors) {
			for (Node gat : gateways) {
				for (int sf : SF) {
					//double val = (Hss.get(i.getName()).get(gat.getName()).get(sf) - beta) + 0.99;
					double val = Math.ceil(Hss.get(i.getName()).get(gat.getName()).get(sf) - beta);
					IloLinearNumExpr zicns = cplex.linearNumExpr();
					// Hss = new HashMap<String, HashMap<String, HashMap<Integer, Double>>>();
					zicns.addTerm(j.get(i).get(gat).get(sf), val);
					//cplex.addLe(z.get(i).get(gat).get(sf), zicns);
					cplex.addEq(z.get(i).get(gat).get(sf), zicns);

				}
			}
		}
		// Zikf geq constraint
		for (Node i : sensors) {
			IloLinearNumExpr zcns = cplex.linearNumExpr();
			for (Node gat : gateways) {
				for (int sf : SF) {
					zcns.addTerm(1, z.get(i).get(gat).get(sf));

				}
			}
			cplex.addGe(zcns, 1);
		}

		if (mode == 1) {
			// Linearization for W//y,gat,j,sf1,sf2
			for (Node i : sensors) {
				for (Node gat : gateways) {
					for (Node x : sensors) {
						for (int sf : SF) {
							for (int sf1 : SF) {
								if (!i.equals(x)) {
									cplex.addLe(w.get(i).get(gat).get(x).get(sf).get(sf1), j.get(i).get(gat).get(sf));
									cplex.addLe(w.get(i).get(gat).get(x).get(sf).get(sf1), y.get(x).get(sf1));
									IloLinearNumExpr lincns = cplex.linearNumExpr();
									lincns.addTerm(1, j.get(i).get(gat).get(sf));
									lincns.addTerm(1, y.get(x).get(sf1));
									IloNumVar n = cplex.numVar(-1, 1);
									lincns.addTerm(1, n);
									cplex.addGe(w.get(i).get(gat).get(x).get(sf).get(sf1), lincns);

								}
							}
						}
					}
				}
			}

			double val = k.size();
			// A2 constraint
			for (Node i : sensors) {
				for (Node j : sensors) {
					for (int sf : SF) {

						IloLinearNumExpr acns = cplex.linearNumExpr();
						for (Node gat : gateways) {
							acns.addTerm(1 * (1 / val), k.get(gat));
							if (!i.equals(j))
								acns.addTerm(-(1 / val) * Cf.get(i.getName() + j.getName()).get(gat.getName()).get(sf),
										w.get(i).get(gat).get(j).get(sf).get(sf));
						}
						IloLinearNumExpr acns2 = cplex.linearNumExpr();
						IloNumVar n = cplex.numVar(1, 1);
						acns2.addTerm(1, n);
						acns2.addTerm(-1, a.get(i).get(j).get(sf));
						cplex.addGe(acns2, acns);

					}
				}
			}
			// A1 constraint
			for (Node i : sensors) {
				for (Node j : sensors) {
					for (int sf : SF) {

						IloLinearNumExpr acns = cplex.linearNumExpr();
						for (Node gat : gateways) {
							acns.addTerm(1 * (1 / val), k.get(gat));
							if (!i.equals(j))
								acns.addTerm(-(1 / val) * Cf.get(i.getName() + j.getName()).get(gat.getName()).get(sf),
										w.get(i).get(gat).get(j).get(sf).get(sf));
						}
						IloLinearNumExpr acns2 = cplex.linearNumExpr();
						IloNumVar n = cplex.numVar(1, 1);
						IloNumVar n1 = cplex.numVar(1, 0.999999999999999999);
						acns2.addTerm(1, n);
						acns2.addTerm(-1, a.get(i).get(j).get(sf));
						acns.addTerm(1, n1);
						cplex.addLe(acns2, acns);

					}
				}
			}
			//B constr
			for (Node s : sensors) {

				for (Node s2 : sensors) {

					for (int sf : SF) {
						IloLinearNumExpr bcns = cplex.linearNumExpr();
						for (int sf1 : SF) {
							if (sf != sf1) {
								for (Node gat : gateways) {

									bcns.addTerm(1 * (1 / val), k.get(gat));
									if (!s.equals(s2)) {
										bcns.addTerm(
												-(1 / val) * If.get(s.getName() + s2.getName()).get(gat.getName())
														.get(sf + "_" + sf1),
												w.get(s).get(gat).get(s2).get(sf).get(sf1));
									}
								}
								IloLinearNumExpr bcns2 = cplex.linearNumExpr();
								IloNumVar n = cplex.numVar(1, 1);
								IloNumVar n1 = cplex.numVar(1, 0.999999999999999999);
								bcns2.addTerm(1, n);
								bcns2.addTerm(-1, b.get(s).get(s2).get(sf).get(sf1));
								bcns.addTerm(1, n1);
								cplex.addLe(bcns2, bcns);

							}
						}

					}

				}

			}
			for (Node s : sensors) {

				for (Node s2 : sensors) {

					for (int sf : SF) {
						IloLinearNumExpr bcns = cplex.linearNumExpr();
						for (int sf1 : SF) {
							if (sf != sf1) {
								for (Node gat : gateways) {

									bcns.addTerm(1 * (1 / val), k.get(gat));
									if (!s.equals(s2)) {
										bcns.addTerm(
												-(1 / val) * If.get(s.getName() + s2.getName()).get(gat.getName())
														.get(sf + "_" + sf1),
												w.get(s).get(gat).get(s2).get(sf).get(sf1));
									}
								}
								IloLinearNumExpr bcns2 = cplex.linearNumExpr();
								IloNumVar n = cplex.numVar(1, 1);
								//IloNumVar n1 = cplex.numVar(1, 0.999999999999999999);
								bcns2.addTerm(1, n);
								bcns2.addTerm(-1, b.get(s).get(s2).get(sf).get(sf1));
								//bcns.addTerm(1, n1);
								cplex.addGe(bcns2, bcns);

							}
						}

					}

				}

			}

		// ultimo vincolo con solo A
		// (ToA[sf-7]/1000.0)-Math.log(threshold)/(2*lambda)
		
		for (Node i : sensors) {
			IloLinearNumExpr acount = cplex.linearNumExpr();
			for (Integer sf : SF) {
				//IloLinearNumExpr acount = cplex.linearNumExpr();
				// IloLinearNumExpr acountf = cplex.linearNumExpr();
				for (Node j : sensors) {
					if (!i.equals(j)) {
						acount.addTerm((ToA[sf - 7] / 1000.0), a.get(i).get(j).get(sf));
					}
				}
				
				// IloNumVar n = cplex.numVar(1, ToA[sf-7]);
				acount.addTerm((ToA[sf - 7] / 1000.0), cplex.numVar(1, 1));
				//cplex.addLe(acount, Math.log(threshold) / (-2 * lambda)); // Fixed
				

			}
			for(Integer sf : SF)
			{
				//IloLinearNumExpr acount = cplex.linearNumExpr();
				for(Integer sf1: SF)
				{
					for(Node j :sensors)
					{
						if(sf!=sf1)
						{
							if(!i.equals(j))
							{
								acount.addTerm((ToA[sf - 7] / 1000.0), b.get(i).get(j).get(sf).get(sf1));
							}
						}
					}
				}
			}
			cplex.addLe(acount, Math.log(threshold) / (-2 * lambda)); // Fixed

		}
		//COPIA
		/*for (Node i : sensors) {
			for (Integer sf : SF) {
				IloLinearNumExpr acount = cplex.linearNumExpr();
				// IloLinearNumExpr acountf = cplex.linearNumExpr();
				for (Node j : sensors) {
					if (!i.equals(j)) {
						acount.addTerm((ToA[sf - 7] / 1000.0), a.get(i).get(j).get(sf));
					}
				}
				
				// IloNumVar n = cplex.numVar(1, ToA[sf-7]);
				acount.addTerm((ToA[sf - 7] / 1000.0), cplex.numVar(1, 1));
				cplex.addLe(acount, Math.log(threshold) / (-2 * lambda)); // Fixed
				

			}

		}*/
	}

	}

	public void solve(String f) throws IloException, IOException {
		// long resolTime = System.currentTimeMillis();
		if (cplex.solve()) {
			BufferedWriter wout = new BufferedWriter(new FileWriter(f));
			BufferedWriter wout7 = new BufferedWriter(new FileWriter("nodisf7.csv"));
			BufferedWriter wout8 = new BufferedWriter(new FileWriter("nodisf8.csv"));
			BufferedWriter wout9 = new BufferedWriter(new FileWriter("nodisf9.csv"));
			BufferedWriter wout10 = new BufferedWriter(new FileWriter("nodisf10.csv"));
			BufferedWriter wout11 = new BufferedWriter(new FileWriter("nodisf11.csv"));
			BufferedWriter wout12 = new BufferedWriter(new FileWriter("nodisf12.csv"));
			BufferedWriter wout2 = new BufferedWriter(new FileWriter("gateways.csv"));
			BufferedWriter wout3 = new BufferedWriter(new FileWriter("gatewaysn.csv"));
			wout.write("Solved \n");
			wout7.write("C00" + "," + "0.0" + "," + "0.0" + "," + "7" + "\n");
			wout8.write("C00" + "," + "0.0" + "," + "0.0" + "," + "8" + "\n");
			wout9.write("C00" + "," + "0.0" + "," + "0.0" + "," + "9" + "\n");
			wout10.write("C00" + "," + "0.0" + "," + "0.0" + "," + "10" + "\n");
			wout11.write("C00" + "," + "0.0" + "," + "0.0" + "," + "11" + "\n");
			wout12.write("C00" + "," + "0.0" + "," + "0.0" + "," + "12" + "\n");
			wout3.write("0.0"+","+"0.0"+"\n");
			for (Node y1 : y.keySet()) {
				for (int sf : y.get(y1).keySet()) {
					if (cplex.getValue(y.get(y1).get(sf)) > eps) {
						wout.write("Node " + y1.getName() + " is allocated with Sf " + sf + "\n");
						System.out.println("Gateway " + y1.getName() + " is allocated with Sf " + sf + "\n");
						switch (sf) {
						case 7:
							wout7.write(y1.getName() + "," + y1.getX() + "," + y1.getY() + "," + sf + "\n");
							break;
						case 8:
							wout8.write(y1.getName() + "," + y1.getX() + "," + y1.getY() + "," + sf + "\n");
							break;
						case 9:
							wout9.write(y1.getName() + "," + y1.getX() + "," + y1.getY() + "," + sf + "\n");
							break;
						case 10:
							wout10.write(y1.getName() + "," + y1.getX() + "," + y1.getY() + "," + sf + "\n");
							break;
						case 11:
							wout11.write(y1.getName() + "," + y1.getX() + "," + y1.getY() + "," + sf + "\n");
							break;
						default:
							wout12.write(y1.getName() + "," + y1.getX() + "," + y1.getY() + "," + sf + "\n");
						}
					}
				}

			}
			wout.write("Gateway usati :\n");
			for (Node gat : k.keySet()) {
				System.out.println("Valore di " + gat.getName() + " " + cplex.getValue(k.get(gat)));
				if (cplex.getValue(k.get(gat)) == 1) {
					wout2.write(gat.getX() + "," + gat.getY() + "\n");
					wout.write("Gat " + gat.getName() + "\n");
				} else {
					wout3.write(gat.getX() + "," + gat.getY() + "\n");
				}

			}
			wout.close();
			wout7.close();
			wout8.close();
			wout9.close();
			wout10.close();
			wout11.close();
			wout12.close();
			wout2.close();
			wout3.close();
			for (Node i : sensors) {
				for (Node j : sensors) {
					for (int sf : SF) {
						if (!i.equals(j)) {
							// System.out.println(a.get(i).get(j).get(sf) + " = " +
							// cplex.getValue(a.get(i).get(j).get(sf)));
						}

					}
				}
			}
			// Map<Node, Map<Node, Map<Integer, Map<Integer, IloNumVar>>>> w1 = new
			// HashMap<Node, Map<Node, Map<Integer, Map<Integer, IloNumVar>>>>();
			for (Node i : sensors) {
				for (Node gat : gateways) {
					for (Node j : sensors) {
						for (int sf : SF) {
							for (int sf1 : SF) {
								if (!i.equals(j)) {
									/*
									 * System.out.println("i" + i.getName() + "j" + j.getName() + "gat " +
									 * gat.getName() + "  "+ sf + " " + sf1 + " "
									 * +cplex.getValue(w.get(i).get(gat).get(j).get(sf).get(sf1)));
									 */
								}
							}
						}
					}
				}

			}
		}
	}

	boolean isServed(Node i) throws UnknownObjectException, IloException {
		for (int sf : y.get(i).keySet())
			if (cplex.getValue(y.get(i).get(sf)) > eps)
				return true;
		return false;
	}

	public void setLP() throws IloException {

		// s = new HashMap<Node, Map<Node, Map<Integer, IloNumVar>>>();

		// Variable Sif
		for (Node cli : sensors) {
			HashMap<Integer, IloNumVar> s1 = new HashMap<Integer, IloNumVar>();

			for (Integer sf : SF) {
				s1.put(sf, cplex.numVar(0, 1));
			}
			s.put(cli, s1);

		}

		// R for the linearization of Sijk and yif
		// new HashMap<Node, Map<Integer, Map<Node, Map<Integer, IloNumVar>>>>();
		for (Node y : sensors) {
			HashMap<Integer, Map<Node, Map<Integer, IloNumVar>>> r1 = new HashMap<Integer, Map<Node, Map<Integer, IloNumVar>>>();
			for (int sf : SF) {
				HashMap<Node, Map<Integer, IloNumVar>> r2 = new HashMap<Node, Map<Integer, IloNumVar>>();
				for (Node j : sensors) {
					if (y.equals(j)) {
						HashMap<Integer, IloNumVar> r3 = new HashMap<Integer, IloNumVar>();
						for (int sf1 : SF) {
							if (sf1 == sf)
								r3.put(sf1, cplex.numVar(0, 1));
						}
						r2.put(j, r3);
					}

				}
				r1.put(sf, r2);
			}
			r.put(y, r1);
		}

		// Settings variables for gateways Xk.
		for (Node gat : gateways)
			k.put(gat, cplex.boolVar("K" + gat.getName()));

		// Setting variable y^sf_i
		for (Node s : sensors) {
			Map<Integer, IloNumVar> ys = new HashMap<Integer, IloNumVar>();

			for (int sf : SF) // Setting variables for clients
				ys.put(sf, cplex.boolVar("y^" + sf + "_" + s.getName()));

			y.put(s, ys);
			// System.out.println(s.getName()+" "+y.get(s));
		}
		// Setting the variable JXkYif
		for (Node s : sensors) {
			Map<Node, Map<Integer, IloNumVar>> j2 = new HashMap<Node, Map<Integer, IloNumVar>>();
			for (Node gat : gateways) {
				Map<Integer, IloNumVar> j1 = new HashMap<Integer, IloNumVar>();
				for (int sf : SF) {
					j1.put(sf, cplex.boolVar("X_" + gat.getName() + "y^" + "sf" + "_" + s.getName() + sf));
				}

				// put in interm
				j2.put(gat, j1);
			}
			// put in j
			j.put(s, j2);
		}

		// Setting the variable Z

		for (Node s : sensors) {
			Map<Node, Map<Integer, IloNumVar>> z2 = new HashMap<Node, Map<Integer, IloNumVar>>();
			for (Node gat : gateways) {
				Map<Integer, IloNumVar> z1 = new HashMap<Integer, IloNumVar>();
				for (int sf : SF) {
					z1.put(sf, cplex.boolVar("Z_y" + s.getName() + " k_" + gat.getName() + " sf_" + sf));
				}
				z2.put(gat, z1);
			}
			z.put(s, z2);

		}

		// Setting the variable W-JXkYifyjf
		for (Node y : sensors) {
			Map<Node, Map<Node, Map<Integer, Map<Integer, IloNumVar>>>> w1 = new HashMap<Node, Map<Node, Map<Integer, Map<Integer, IloNumVar>>>>();
			for (Node gat : gateways) {
				Map<Node, Map<Integer, Map<Integer, IloNumVar>>> w2 = new HashMap<Node, Map<Integer, Map<Integer, IloNumVar>>>();
				for (Node j : sensors) {
					Map<Integer, Map<Integer, IloNumVar>> w3 = new HashMap<Integer, Map<Integer, IloNumVar>>();
					if (!(y.equals(j))) {
						for (int sf1 : SF) {
							Map<Integer, IloNumVar> w4 = new HashMap<Integer, IloNumVar>();
							for (int sf2 : SF) {
								w4.put(sf2, cplex.boolVar(
										"X_" + gat.getName() + y.getName() + "y^" + sf1 + j.getName() + "j^ " + sf2));
							}
							w3.put(sf1, w4);
						}
					}
					w2.put(j, w3);
				}
				w1.put(gat, w2);

			}
			w.put(y, w1);
		}
		// setting variable A
		for (Node s : sensors) {
			Map<Node, Map<Integer, IloNumVar>> a1 = new HashMap<Node, Map<Integer, IloNumVar>>();
			for (Node s2 : sensors) {
				Map<Integer, IloNumVar> a2 = new HashMap<Integer, IloNumVar>();
				for (int sf : SF) {
					a2.put(sf, cplex.boolVar("i_" + s.getName() + "j_" + s2.getName() + "sf" + sf));
				}
				a1.put(s2, a2);
			}
			a.put(s, a1);
		}

		// setting variable B
		for (Node s : sensors) {
			Map<Node, Map<Integer, Map<Integer, IloNumVar>>> b1 = new HashMap<Node, Map<Integer, Map<Integer, IloNumVar>>>();
			for (Node s2 : sensors) {
				Map<Integer, Map<Integer, IloNumVar>> b2 = new HashMap<Integer, Map<Integer, IloNumVar>>();
				for (int sf : SF) {
					Map<Integer, IloNumVar> b3 = new HashMap<Integer, IloNumVar>();
					for (int sf1 : SF) {
						if (sf != sf1) {
							b3.put(sf1,
									cplex.boolVar("i_" + s.getName() + "j_" + s2.getName() + "sf" + sf + "sf1" + sf1));
						}
					}
					b2.put(sf, b3);
				}
				b1.put(s2, b2);
			}
			b.put(s, b1);
		}
	}
}
