///4lm4n4ch,envagyokalevente@gmail.com


import java.util.*;


import game.racetrack.*;
import game.racetrack.utils.*;

import static java.lang.Math.*;

/**
 * Az Agent osztály, ami az A* algoritmust használja a legrövidebb út meghatározására, figyelembe véve a pályán található akadályokat és érméket.
 */
public class Agent extends RaceTrackPlayer {
    // Célvonal és kezdőpont koordinátái
    int targetI, targetJ,startI,startJ;
    int asd=0;

    PriorityQueue<Node> openSet; //pontok halmaza, amiből mindig a legkisebb fCostal rendelkezőt nézzük meg
    HashSet<Node> closedSet; //pontok, amiket már megnéztünk
    Stack<Node> path = new Stack<>(); //a legjobb utat tartalmazó Stack
    HashMap<Node , Node> openSetMap = new HashMap<>();


    /**
     * Konstruktor, amelyel létrehozzuk az osztályt.
     * @param state A játékos aktuális állapota.
     * @param random egy random szám, amit használhat a kód.
     * @param track A versenypálya.
     * @param coins Az érmék tömbbe rendezve.
     * @param color A játékos szine.
     */
    public Agent(PlayerState state, Random random, int[][] track, Coin[] coins, int color) {
        super(state, random, track, coins, color);

        // Inicializálás a konstruktorban
        openSet = new PriorityQueue<>(Comparator.comparingLong(n -> n.fCost));
        closedSet = new HashSet<>();

        //System.out.println(isWall(state.i,state.j) +"illetve fale: "+ isWall(0,0));
        //startpont inicializálás
        startI = state.i;
        startJ = state.j;
        //System.out.println("kezdés: " + startI +"]["+ startJ);
        //cél meghatározás
        setFinishLine();
        aalgorithm();
        Node targetNode = null;
        for (Node node : closedSet) {
            //System.out.println(targetI +"][" + targetJ +":  " + node.toString());
            if (isFinished(node)) {
                //System.out.println("magtalálta  a végét: " + node.toString());
                targetNode = node;
                break;
            }
        }

        if (targetNode != null) {
            Node current = targetNode;
            while (current != null) {
                //System.out.println("A list: " + current.toString());
                path.push(current);
                current = current.parent;
            }
            // Töröljük a startpontot, mivel már ott vagyunk
            if (!path.isEmpty()) {
                path.pop();
            }
        }
    }
    /**
     * Visszaadja a következő lépés irányát a játékos számára.
     * @param remainingTime A rendelkezésre álló idő.
     * @return Az új irány, amerre a játékosnak mozognia kell.
     */
    @Override
    public Direction getDirection(long remainingTime) {
        //if(true) return null;
        Node n = path.pop();
        //System.out.println("itt az új irány: " + n.toString());
        if(true) return new Direction(n.i-n.parent.i-state.vi,n.j-n.parent.j-state.vj);
        return pathFinder();
    }
    /**
     * Az A* algoritmus implementációja, optimális útvonal megtalálása.
     */
    public void aalgorithm(){
        //start Node openhez adás
        openSet.add(new Node(startI, startJ, null));

        int h=-10;
        while(!openSet.isEmpty()){ //addig fusson, amíg van az opensetben.
            //current Node az openSet legalacsonyabb f-costja
            Node currentNode = openSet.poll();
            //System.out.println("Jelenlegi csomópont: (" + currentNode.i + ", " + currentNode.j + ") fCost: " + currentNode.fCost);

            closedSet.add(currentNode); //aktuális Node closSet be rakása
            openSet.remove(currentNode);

            //letesztelni, hogy megoldás e az aktuális pont
            if(isFinished(currentNode));
            // A currentNode összes szomszédjának feldolgozása
            List<Node> neighbors = listNeighbors(currentNode);
            for (Node neighbor : neighbors) {
                //System.out.println(neighbor.toString() +"\nSzülője  " +neighbor.parent.toString()+"\n\n");
                if (closedSet.contains(neighbor)) {
                    continue; // A szomszéd már szerepel a zárt halmazban
                }

                if (!openSet.contains(neighbor)) {
                    openSet.add(neighbor); // Ha még nem szerepel a nyílt halmazban, akkor hozzáadjuk
                } else {
                    // Ha már szerepel az openSet-ben, frissítjük a költségeket, ha olcsóbb útvonalat találtunk

                    Node existing = openSet.stream().filter(n -> n.equals(neighbor)).findFirst().orElse(null);
                    if (existing != null && neighbor.fCost < existing.fCost) {
                        existing.parent = currentNode;
                        existing.gCost = neighbor.gCost;
                        existing.fCost = neighbor.fCost;
                        existing.hCost = neighbor.hCost;
                        existing.vi = neighbor.vi;
                        existing.vj = neighbor.vj;
                    }
                }
            }
        }
    }
    /**
     * Kiszámítja a heurisztikus költséget a céltól való távolság alapján.
     * Figyelembe veszi az érméket is, amennyiben azok az adott ponton találhatóak.
     * @param toI A cél csomópont sor indexe.
     * @param toJ A cél csomópont oszlop indexe.
     * @return A heurisztika költsége.
     */
    private long calculateHeuristic(int toI, int toJ){
        int di = abs(toI - targetI);
        int dj = abs(toJ - targetJ);
        double cost;
        int diagonalSteps = min(di, dj);
        int straightSteps = di + dj - 2 * diagonalSteps; //a keresztbe lépés után,m ebből megadjuk az egyenes lépéseket, majd ebből kiszámoljuk a heurisztikát
        cost= 16 * diagonalSteps + 10 * straightSteps;
        return (long) cost;
    }
    /**
     * Beállítja a célvonalat a pályán. A célvonal koordinátáit a 'targetI' és 'targetJ' változókban tárolja.
     */
    public void setFinishLine() {
        for (int i = 0; i < super.track.length; i++) {
            for (int j = 0; j < super.track[i].length; j++) {
                if ((super.track[i][j] & RaceTrackGame.FINISH) == RaceTrackGame.FINISH) {
                    targetI = i;
                    targetJ = j;
                    return;
                }
            }
        }
       // System.err.println("Nem találja a célt!!"); //cél beállítása, ha nem találja
        targetI = 0;
        targetJ = super.track[0].length - 1;
    }
    /**
     * Ellenőrzi, hogy a megadott csomópont a célvonalon található-e.
     * @param node A vizsgálandó Node.
     * @return Igaz, ha célba ért van, egyébként hamis.
     */
    public boolean isFinished(Node node) {
        return (super.track[node.i][node.j] & RaceTrackGame.FINISH) == RaceTrackGame.FINISH;
    }

    /**
     * Ellenőrzi, hogy a megadott koordinátákon fal található-e.
     * @param i A sor indexe.
     * @param j Az oszlop indexe.
     * @return Igaz, ha fal van a megadott helyen, egyébként hamis.
     */
    private boolean isWall(int i, int j) {
        return (track[i][j] & RaceTrackGame.WALL) == RaceTrackGame.WALL;
    }
    /**
     * Ellenőrzi, hogy a megadott koordináták érvényesek-e a pályán belül.
     * @param i A sor indexe.
     * @param j Az oszlop indexe.
     * @return Igaz, ha a koordináták érvényesek a pályán, egyébként hamis.
     */
    private boolean isValidPosition(int i, int j) {
        return i >= 0 && i < track.length && j >= 0 && j < track[0].length;
    }
    /**
     * Kiszámítja a G költséget egy adott csomópontból egy másikba való mozgáshoz.
     * Figyelembe veszi a sebességet és az érmék jelenlétét.
     * @param from A kiinduló csomópont.
     * @param toI A cél csomópont sor indexe.
     * @param toJ A cél csomópont oszlop indexe.
     * @return A G költség értéke.
     */
    private long calculateGCost(Node from, int toI, int toJ) {
        if (from.parent == null) {
            from.velocity = 0; // Kezdeti sebesség
            return 0;
        }
        // Sebesség számítása
        if (isSameDirection(from, toI, toJ)) {
            from.velocity++;
        } else {
            from.velocity = 1;
        }
        double cost=from.gCost;
        //sebesség beleszámítása
        double speedFactor = calculateSpeedFactor(from.velocity);
        // Érmék értékének figyelembe vétele
        if (isCoinAt(toI, toJ)) {
            cost *= 0.7;
        }
        if (from.i == from.parent.i || from.j == from.parent.j) {
            cost += 10 - speedFactor;
        } else {
            cost += 14 - speedFactor;
        }

        return (long) cost;
    }
    /**
     * Ellenőrzi, hogy van-e érme a megadott koordinátákon.
     * @param i A sor indexe.
     * @param j Az oszlop indexe.
     * @return Igaz, ha érme található a megadott helyen, egyébként hamis.
     */
    private boolean isCoinAt(int i, int j) {
        return (track[i][j] & RaceTrackGame.COIN) == RaceTrackGame.COIN;
    }
    /**
     * Kiszámítja a sebesség szorzót, ami a G költségbe kerül beleszámításra.
     * @param velocity Az aktuális sebesség értéke.
     * @return A sebesség szorzó értéke.
     */
    private double calculateSpeedFactor(int velocity) {
        return min(velocity, 5);
    }
    /**
     * Ellenőrzi, hogy a csomópont ugyanabban az irányban halad-e, mint az előző lépésben.
     * @param node Az aktuális csomópont.
     * @param toI A cél csomópont sor indexe.
     * @param toJ A cél csomópont oszlop indexe.
     * @return Igaz, ha az irány megegyezik, egyébként hamis.
     */
    private boolean isSameDirection(Node node, int toI, int toJ) {
        int dirI = toI - node.i;
        int dirJ = toJ - node.j;
        int prevDirI = node.i - node.parent.i;
        int prevDirJ = node.j - node.parent.j;
        return (dirI == prevDirI && dirJ == prevDirJ);
    }
    /**
     * Lekérdezi a szomszédos csomópontok listáját az aktuális csomóponthoz képest.
     * Kizárja azokat a szomszédokat, amelyek falak vagy érvénytelen pozíciók, vagy az előző lépéssel ellentétes irányba haladnak.
     * @param currentNode Az aktuális csomópont.
     * @return A szomszédos csomópontok listája.
     */
    public List<Node> listNeighbors(Node currentNode) {
        List<Node> neighbors = new ArrayList<>();
        for (Direction dir : RaceTrackGame.DIRECTIONS) {

            int ci = currentNode.vi + dir.i; //szomszéd irányának koordinátái
            int cj = currentNode.vj + dir.j; //
            Direction directionWithVelocity = setDirectionbyVelocity(dir);
            boolean isGood=true;
            if(0>currentNode.i || currentNode.i+ci>=track.length || 0 > currentNode.j || currentNode.j>= track[currentNode.i].length){
                //System.out.println("pályám kívül");
                isGood= false; //ha pályán kívül van
            }
            if(ci == 0 && cj ==0) { //its good
                //System.out.println("ugyanaz");
                isGood=false;
            }
            int maxCiAbs = max(abs(ci),abs(cj));
            int minCiAbs = min(abs(ci),abs(cj));
            //végigmenni, hogy van e fal a kettő között

            int m,n;
            if(ci<0) m=-1; else m=1;
            if(cj<0) n=-1; else n=1;

            Node futoN = new Node(
                    currentNode.i, currentNode.j,currentNode.parent); //parentet nem kell használni
            Node targetN = new Node(currentNode.i +ci,currentNode.j+cj,null);

            while(!(futoN.i == targetN.i && futoN.j == targetN.j)){



                //System.out.println(futoN.i +"=="+ targetN.i +"&&"+ futoN.j+" == "+targetN.j);
                if(m*(targetN.i-futoN.i) > n*(targetN.j-futoN.j))  {
                   // System.out.println("i nagyobb");
                    futoN.i+=m;
                }
                else if((targetN.i-futoN.i)*m<n*(targetN.j-futoN.j)) {
                   // System.out.println("j nagyobb");
                    futoN.j+=n;
                }
                else {
                  // System.out.println("ugyanannyi");
                    futoN.i+=m;
                    futoN.j+=n;
                }


                if(!RaceTrackGame.isNotWall(futoN.i,futoN.j,track)) {

                    //System.out.println("its a wall: " +futoN.i + "][" + futoN.j);
                    isGood=false;
                    break;
                }

            }
            futoN = new Node(
                    currentNode.i, currentNode.j,currentNode.parent); //parentet nem kell használni
            while(!(futoN.i == targetN.i && futoN.j == targetN.j)){
                //System.out.println(futoN.i +"=="+ targetN.i +"&&"+ futoN.j+" == "+targetN.j);
                if(targetN.i!=futoN.i)  futoN.i+=m;
                if(targetN.j!=futoN.j) futoN.j+=n;

                if(!RaceTrackGame.isNotWall(futoN.i,futoN.j,track)) {

                    //System.out.println("its a wall: " +futoN.i + "][" + futoN.j);
                    isGood=false;
                    break;
                }

            }



            if (!isGood) continue; //maybe wrong place


            Node neighbour = new Node(currentNode.i +ci,currentNode.j+cj,ci,cj,currentNode);
            if(neighbour.parent!=null) neighbour.gCost  = neighbour.parent.gCost+12; //test if null
            else neighbour.gCost=10;
            neighbour.hCost = calculateHeuristic(currentNode.i +ci,currentNode.j+cj);

            if(isCoinAt(neighbour.i, neighbour.j)) {
                double cost= neighbour.gCost - getCoinValueAt(neighbour.i, neighbour.j);
                neighbour.gCost = (long)cost;
            }
            neighbour.fCost = neighbour.gCost+ neighbour.hCost;

            neighbors.add(neighbour);

        }
        //System.out.println("Mérete: "+neighbors.size());

        return neighbors;
    }
    public int getCoinValueAt(int i, int j) {
        for (Coin coin : coins) {
            if (coin.i == i && coin.j == j) {
                return coin.value; // Feltételezve, hogy van egy 'value' mező a Coin osztályban
            }
        }
        return 0; // Ha nincs Coin ezen a pozíción, akkor 0-t adunk vissza
    }

    class Node {

        int i; //node sor(!) poziciója
        int j; //node oszlop(!) pozíciója

        int vi;
        int vj;

        long gCost;  //távolság a kezdőponttól
        long hCost;  //távolság a végétől
        long fCost;  //f=g+h

        int velocity=0;  //outdated
        Direction direction; //outdated

        Node parent; //szűlő

        public Node(int i, int j, Node parent, long gCost, long hCost) { //konstruktor
            this.i = i;
            this.j = j;
            this.parent = parent;
            this.gCost = gCost;
            this.hCost = hCost;
            this.fCost = gCost + hCost;
        }
        public Node(int i, int j, Node parent) { //konstruktor heurisztika nélküli létrehozás
            this.i = i;
            this.j = j;
            this.parent = parent;
            //G cost kiszámítás
            this.gCost = 0;
            //heurisztika kiszámítás
            this.hCost = 0;

            this.fCost = gCost+hCost;
        }

        public Node(int i, int j, int vi, int vj, Node parent) {
            this.i = i;
            this.j = j;
            this.vi = vi;
            this.vj = vj;
            this.parent = parent;
        }
        public void setDirection(){
            if (parent != null) {
                int diffX =(int) signum(this.i - parent.i);
                int diffY =(int) signum(this.j - parent.j);

                // A Direction objektum létrehozása az i és j irányértékekkel
                // Az irányvektor értékei -1, 0, vagy 1 lehetnek, a különbség előjelétől függően
                this.direction = new Direction(diffX, diffY);
            } else {
                // A kezdő Node esetén az irány értéke 0 lehet mindkét tengely mentén
                this.direction = new Direction(0, 0);
            }
        }

        @Override
        public boolean equals(Object o) {
            //setDirection();
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Node node = (Node) o;
            //System.out.println(i + "eq " + node.i + "and " + j+ "eq " + node.j);
            return i == node.i && j == node.j;
        }

        @Override
        public int hashCode() {
            //setDirection();
            return Objects.hash(i, j, vi,vj);
        }

        @Override
        public String toString() {
            return "Node{" +
                    "i=" + i +
                    ", j=" + j +
                    ", vi=" + vi +
                    ", vj=" + vj +
                    ", gCost=" + gCost +
                    ", hCost=" + hCost +
                    ", fCost=" + fCost +
                    '}';
        }
    }

    /**
     * Meghatározza az ideális irányt.
     * @return A következő lépés iránya.
     */
    public Direction pathFinder(){

        Direction velocity = new Direction(state.vi,state.vj); //aktuális sebesség
        Direction direction = setDirectionbyVelocity(velocity);
        //kiszámoljuk mennyi pontot ugrottunk át, és annyit kiveszünk belőle

        //if(true) return null;
        Node nextNode = path.peek(); //Kövi elem
        int di = nextNode.i - state.i;
        int dj = nextNode.j - state.j;
        //System.out.println("state v: " + state.i +"]["+ state.j +" nextNode: " +nextNode.i +"]["+ nextNode.j);
        //végigézi a path-ot, addig, amíg egyenesen meg, megszámolja

        //ha nagyobb a táv, mint a sebesség, akkor gyorsítjon,

        //ha megvan a direction, amennyir megy előre vegyen ki a path-ból




        //System.out.println(direction.i +" és " +direction.j+ "velocity " +state.vi + state.vj);

        //mivel ha átlósan megy, a sebességnek ugyanannyinak kell lennie, ezt letesztelni
        Stack<Node> cpath = new Stack<>();
        cpath.addAll(path);
        //meddig megy egyenesen
        int futof=0;
        //if(true)return truckSimulator(di,dj);
        Node curent=new Node(
                state.i,
                state.j,
                null);

        while (!cpath.empty()){
            Node next = cpath.pop();//kivesszük a következőt a sorból
            Direction nextDirection = directionByNodes(curent,next);

            if(nextDirection.i == direction.i && nextDirection.j == direction.j) {//ha egy irányba van a sebességgel
                //System.out.println("Egyenesen tovább");
                futof++; //meddig megy egyenesen futof be számolja

            } else break;//todo
            curent=next;
        }
        //System.out.println(di + "][" +dj + "\n next in line " + path.peek().toString());
        int d = max(abs( state.i - curent.i),  abs(state.j-curent.j));//todo
        int v = max(abs(state.vi),abs(state.vj));
        if(futof==0) { //patika
            //if(path.peek().i == state.i && path.peek().j == state.j) path.pop();
            path.pop();
            return new Direction(di-state.vi,dj-state.vj);
        }

        if(f(d,v+1)) { //patkia
            for (int i=0;i<=v;i++) path.pop();
            //System.out.println("Speed Up");
            return new Direction(di,dj);
        }
        if(f(d,v) ) {
            for (int i=0;i<v;i++)
                path.pop();
            //System.out.println("Speed stay");
            return new Direction(0,0); //todo elseif
        }
        for (int i=0;i<v-1;i++) path.pop();
        //System.out.println("Speed down");
        return new Direction(-state.vi,-state.vj);
    }
    /**
     * Megnézi, hogy le tudna e lassulni az adott sebességről, adott idő alatt
     * @return lassulás hossza.
     */
    public boolean f(int d, int v){
        //System.out.println("D éréke: " + d + ", és v értéke: " + v);
        int k=0;
        for (int i=0;i<=v;i++){
            k+=i;
        }
        return d >= k;
    }
    /**
     * Kiszámítja az irányt két csomópont között.
     * @param a Az első csomópont.
     * @param b A második csomópont.
     * @return Az irány a két csomópont között.
     */
    public Direction directionByNodes(Node a,Node b){
        return new Direction(b.i-a.i,b.j-a.j);
    }

    /**
     * Meghatározza a sebesség alapján az irányt.
     * @param v A jelenlegi sebesség.
     * @return Az irány, amely a sebesség alapján van meghatározva.
     */
    public Direction setDirectionbyVelocity(Direction v){
        if(v.i==0 && v.j==0) return new Direction(0,0);
        else if(v.i==0) return new Direction(0, v.j / abs(v.j)); //ha i nulla, akkor
        else if(v.j==0)return new Direction(v.i / abs(v.i),0); //ha j nulla, akkor
        return new Direction( v.i / abs(v.i), v.j / abs(v.j)); //ha nem nulla, akkor
    }
}