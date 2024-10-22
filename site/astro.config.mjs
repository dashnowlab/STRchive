import { defineConfig } from "astro/config";

// https://astro.build/config
export default defineConfig({
  site: "https://strchive.org",
  /** https://github.com/withastro/astro/issues/4190 */
  trailingSlash: "never",
});
