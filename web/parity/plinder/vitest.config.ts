import { defineConfig } from 'vitest/config';
export default defineConfig({test:{include:['parity/plinder/audit.run.ts'],environment:'node',fileParallelism:false}});
