package rpf.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;


public class FindMathcedReads
    {
    private File bamFile=null;
    private SAMFileReader inputSam=null;
    private boolean containsbamRecord=false;//false : the alignment of the returned SAMRecords need only overlap the interval of interest. 
    public FindMathcedReads()
        {

        }

    public void open(File bamFile)
        {
        close();
        this.bamFile=bamFile;
        this.inputSam=new SAMFileReader(this.bamFile);
        this.inputSam.setValidationStringency(ValidationStringency.SILENT);
        }

    public void close()
        {
        if(inputSam!=null)
            {
            inputSam.close();
            }
        this.inputSam=null;
        this.bamFile=null;
        }

    private int scan(String chromosome,int start,int end) throws IOException
        {
        int nCount=0;
        CloseableIterator<SAMRecord> iter=null;
        try {
            iter=this.inputSam.query(chromosome, start, end, this.containsbamRecord);
            while(iter.hasNext())
                {
                SAMRecord rec=iter.next();
                ++nCount;
                }
            return nCount;
            }
        catch (Exception e) {
            throw new IOException(e);
            }
        finally
            {
            if(iter!=null) iter.close();
            }
        }

    public void run(BufferedReader in) throws IOException
        {
        String line;
        while((line=in.readLine())!=null)
            {
            if(line.isEmpty()) continue;
            String tokens[]=line.split("[\t]");
            String chrom=tokens[0];
            int start=Integer.parseInt(tokens[1]);
            int end=Integer.parseInt(tokens[2]);
            int count=scan(chrom,start,end);
            System.err.println(chrom+"\t"+start+"\t"+end+"\t"+count);
            }
        }

    // java -cp sam-1.16.jar:picard-1.16.jar:. BioStar3414  -bam my.sorted.bam
    public static void main(String[] args)
        {
        File bamFile=null;
        FindMathcedReads app;
        try
            {
            app=new FindMathcedReads();
            int optind=0;
            while(optind<args.length)
                {
                if(args[optind].equals("-h"))
                    {
                    return;
                    }
                else if(args[optind].equals("-bam"))
                    {
                    bamFile=new File(args[++optind]);
                    }
                else if(args[optind].equals("--"))
                    {
                    optind++;
                    break;
                    }
                else if(args[optind].startsWith("-"))
                    {
                    System.err.println("Unnown option: "+args[optind]);
                    return;
                    }
                else
                    {
                    break;
                    }
                ++optind;
                }
            if(bamFile==null)
                {
                System.err.println("undefined BamFile");
                return;
                }
            app.open(bamFile);
            if(optind==args.length)
                {
                app.run(new BufferedReader(new InputStreamReader(System.in)));
                }
            else
                {
                while(optind< args.length)
                    {
                    String inputName=args[optind++];
                    BufferedReader in=new BufferedReader(new FileReader(inputName));
                    app.run(in);
                    in.close();
                    }
                }
            app.close();
            }
        catch(Exception err)
            {
            err.printStackTrace();
            }
        }
    }