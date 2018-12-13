import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

public class ConcurrentThreads_opt {

    private static int NUM_WORKER = Integer.parseInt(System.getProperty("worker","100"));

    private static int ITERATION = Integer.parseInt(System.getProperty("iter","2000"));

    private static int CORE_POOL_SIZE = Integer.parseInt(System.getProperty("pool","1"));

    private static boolean SYNC = Boolean.parseBoolean(System.getProperty("sync","false"));

    private static boolean SIZED = Boolean.parseBoolean(System.getProperty("sized", "false"));

    private static boolean PI = Boolean.parseBoolean(System.getProperty("pi", "true"));

    private static long[] preFilledArray;

    public static void main(final String[] args)
            throws Exception {
        System.out.println("worker, iter, pool, sync, sized, pi, totalTime");

        if (CORE_POOL_SIZE > 0) {
            doPool(CORE_POOL_SIZE);
        } else {
            doPool(1);
            doPool(8);
            doPool(16);
        }
        System.exit(0);
    }

    private static void doPool(int poolSize) throws InterruptedException, java.util.concurrent.ExecutionException {
        final ThreadPoolExecutor service = new ThreadPoolExecutor(poolSize, poolSize, 0L,
                TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
        service.prestartAllCoreThreads();

        preFilledArray = new long[ITERATION];
        final Random r = new Random();
        for (int i = 0; i < ITERATION; i++) {
            preFilledArray[i] = r.nextLong();
        }

        long totalTime = 0;
        for (int i = 0; i < 100; i++) {
            totalTime += doActualWork(service);
        }

        System.out.println(String.format("%s, %s, %s, %s, %s, %s, %s",
                NUM_WORKER, ITERATION, poolSize, SYNC, SIZED, PI, totalTime));

        service.shutdownNow();
    }

    private static long doActualWork(final ThreadPoolExecutor service)
            throws InterruptedException, java.util.concurrent.ExecutionException {

        int remain_num_worker = NUM_WORKER;
        int num_worker_limit = java.lang.Math.min(
            Runtime.getRuntime().availableProcessors(),
            service.getPoolSize());
        int cur_num_worker = java.lang.Math.min(remain_num_worker, num_worker_limit);
        long total = 0;

        do {
            final List<Future<Long>> ff = new ArrayList<>();

            for (int i = 0; i < cur_num_worker; i++) {
                final Future<Long> f;
                if (PI) {
                    f = service.submit(new Pi());
                }
                else {
                    f = service.submit(new RunRun());
                }

                ff.add(f);
            }

            for (final Future<Long> fff : ff) {
                total += fff.get();
            }
            remain_num_worker -= cur_num_worker;
            cur_num_worker = java.lang.Math.min(remain_num_worker, num_worker_limit);
        } while (remain_num_worker != 0);

        return total;
    }

    private static class Pi
            implements Callable<Long> {
        @Override
        public Long call()
                throws Exception {

            BigDecimal d = new BigDecimal(0);
            final long start = System.currentTimeMillis();
            for (int i = 0; i < ITERATION; i++) {
                if (SYNC) {
                    synchronized (this) {
                        final BigDecimal next = new BigDecimal(4);
                        final BigDecimal factor = next.divide(new BigDecimal(i * 2 + 1), 1000, RoundingMode.HALF_UP);
                        if (i % 2 == 0) {
                            d = d.add(factor);
                        }
                        else {
                            d = d.subtract(factor);
                        }
                    }
                }
                else {
                    final BigDecimal next = new BigDecimal(4);
                    final BigDecimal factor = next.divide(new BigDecimal(i * 2 + 1), 1000, RoundingMode.HALF_UP);
                    if (i % 2 == 0) {
                        d = d.add(factor);
                    }
                    else {
                        d = d.subtract(factor);
                    }
                }

            }
            return System.currentTimeMillis() - start;
        }
    }

    private static class RunRun
            implements Callable<Long> {

        @Override
        public Long call()
                throws Exception {

            final Map<Long, Long> poor;
            if (SIZED) {
                poor = SYNC ? new PoorMap<>(ITERATION) : new DefaultMap<>(ITERATION);
            }
            else {
                poor = SYNC ? new PoorMap<>() : new DefaultMap<>();
            }

            final long start = System.currentTimeMillis();
            for (int i = 0; i < ITERATION; i++) {
                long p = preFilledArray[i];
                poor.put(p, p);
            }
            return System.currentTimeMillis() - start;
        }
    }

    private static class PoorMap<K, V>
            extends HashMap<K, V> {

        public PoorMap(final int initialCapacity) {
            super(initialCapacity);
        }

        public PoorMap() {
        }

        @Override
        public V put(final K key, final V value) {
            synchronized (this) {
                return super.put(key, value);
            }
        }
    }

    private static class DefaultMap<K, V>
            extends HashMap<K, V> {

        public DefaultMap(final int initialCapacity) {
            super(initialCapacity);
        }

        public DefaultMap() {
        }

        @Override
        public V put(final K key, final V value) {
            return super.put(key, value);
        }
    }

}
