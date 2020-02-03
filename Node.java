import java.util.Set;

import ilog.concert.IloNumVar;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class Node { 
	private double xcoord;
	private double ycoord;
	
	private double ptx;
	private String name;
	//private HashMap<String, Double> Hs = new HashMap<String, Double>(); //mod by me
	//private ArrayList<Double> Hs = new ArrayList <Double>(); //Array list of H calculated on a possible gateway for each node.
	//problema, non posso prendere dal set in ordine gli HS... pensare ad una soluzione alternativa.
	
	public Node(){}
	
	public Node(double x, double y, double p, String n){
		xcoord = x;
		ycoord = y;
		ptx = p;
		name = n;
		
	}
	
	public double getX(){
		return xcoord;
	}
	
	public double getY(){
		return ycoord;
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String n) {
		name = n;
	}
	
	public double getDistance(Node i) {
		return Math.sqrt(Math.pow(getX()-i.getX(), 2.0)+Math.pow(getY()-i.getY(), 2.0));
	}
	
	public String toString(){
		return name +" ("+ xcoord +","+ ycoord +")";
	}
	
	// TX MODEL
	public static double invDb(double x) {
		return Math.pow(10, x/10);
	}
	
	public static double db(double x) {
		return 10*Math.log10(x);
	}
	
	
	// OKUMURA HATA MODEL
	public double getReceivedPower(Node g) {
		return db(invDb(ptx) * invDb(pathLoss(g)));
	} 
	
	public double pathLoss(Node g) {
		double antennaGain = 6;
		double hb = 15;
		double d = getDistance(g)/1000.0;
		return -(69.55+26.16*Math.log10(868)-13.82*Math.log10(hb)+(44.9-6.55*Math.log10(hb))*Math.log10(d)-2*(Math.log10(868/28)*Math.log10(868/28))-5.4)+antennaGain;
	}
	
	public double getH(Node g, double qf) {
		double N = -174+6+10*Math.log10(125e3); // #dBm
		return Math.exp(-invDb(N+qf)/(invDb(getReceivedPower(g))));
	}
	
	/*public HashMap<String, Double> getArrH(Set<Node> gs, double qf) {
		//double N = -174+6+10*Math.log10(125e3); // #dBm
		for (Node g : gs)
			Hs.put(g.name, this.getH(g, qf));
		
		return Hs;
	}*/
	public double getMaxH(Set<Node> gs, double qf) {
		double N = -174+6+10*Math.log10(125e3); // #dBm
		double res = 0.0;
		for (Node g : gs)
			if (getH(g, qf) > res)
				res = getH(g,qf);
		
		return res;
	}
	
	// SF
	public int getMinSF(Set<Node> gs, Double[] qf, double beta) {
		int minSF = 12;
		for (Node g : gs)
			for (int s=7; s<13; s++)
				if (getH(g, qf[s-7]) >= beta && s < minSF)
					minSF = s;
		return minSF;
	}
}
