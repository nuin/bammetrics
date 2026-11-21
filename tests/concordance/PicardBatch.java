import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

/**
 * Runs many Picard command lines in one JVM, for the seed sweeps of
 * test_downsample_sam.py: each line of the file named first is a Picard
 * class of picard.sam and its arguments, separated by tabs. Prints the
 * exit code of each line, one per line, to standard output; Picard's own
 * log goes to standard error. Run as
 * `java -cp picard.jar tests/concordance/PicardBatch.java LINES`, under
 * the Picard queue lock (picard_lock.held).
 */
public class PicardBatch {
    public static void main(String[] args) throws Exception {
        List<String> lines = Files.readAllLines(Path.of(args[0]));
        for (String line : lines) {
            if (line.isEmpty()) {
                continue;
            }
            String[] fields = line.split("\t");
            String[] rest = java.util.Arrays.copyOfRange(fields, 1, fields.length);
            Class<?> tool = Class.forName("picard.sam." + fields[0]);
            picard.cmdline.CommandLineProgram program =
                (picard.cmdline.CommandLineProgram) tool.getDeclaredConstructor().newInstance();
            int code;
            try {
                code = program.instanceMain(rest);
            } catch (Throwable error) {
                System.err.println(error);
                code = 1;
            }
            System.out.println(code);
            System.out.flush();
        }
    }
}
