import { defineConfig } from "astro/config";
import svgr from "vite-plugin-svgr";
import mdx from "@astrojs/mdx";
import react from "@astrojs/react";

// https://astro.build/config
export default defineConfig({
  site: "https://strchive.org",
  integrations: [mdx(), react({})],
  /** https://github.com/withastro/astro/issues/4190 */
  trailingSlash: "never",
  vite: {
    plugins: [
      svgr({
        svgrOptions: {
          /** https://github.com/gregberge/svgr/discussions/770 */
          expandProps: "start",
          svgProps: {
            className: `{props.className ? props.className + " icon" : "icon"}`,
            "aria-hidden": "true",
          },
        },
      }),
    ],
  },
});
