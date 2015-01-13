package fCLIP;

import java.util.PriorityQueue;
import java.util.Date;

public class Test
{
	    public static void main(String[] args)
    {
        PriorityQueue<Integer> queue = 
            new PriorityQueue<Integer>(2);
       
        queue.add(-11);
        queue.add(4);
        queue.add(3);
        queue.add(3);
        queue.add(3);
        queue.add(3);
        queue.add(3);
        queue.add(3);
        queue.add(30);
        System.out.println(new Date().getTime());
        while (queue.size() > 4)
        {
        	
            System.out.println(queue.remove());
        }
        System.out.println(queue);
    }
}


