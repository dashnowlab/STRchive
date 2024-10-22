import { defineConfig } from "astro/config";
import mdx from "@astrojs/mdx";

// https://astro.build/config
export default defineConfig({
  site: "https://strchive.org",
  integrations: [mdx()],
  /** https://github.com/withastro/astro/issues/4190 */
  trailingSlash: "never",
});
