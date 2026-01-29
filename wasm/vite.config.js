import { defineConfig } from 'vite';
import wasm from 'vite-plugin-wasm';

// The base path is trick, for local test by run through command in `.githut/workflows/deploy.yaml`
// the base path needs to be `/` because index.html in dist start at path /
// While at GitHub page deployment, the url path to the assets becomes `https://rs4rse.github.io/vizmat/assets/`.
// I need the base path be `/vizmat` to correctly point to the assets.
export default defineConfig({
  plugins: [wasm()],
  base: '/vizmat',
  optimizeDeps: {
    exclude: ['vizmat']
  }
});
