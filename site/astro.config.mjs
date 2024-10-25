import { defineConfig } from "astro/config";
import react from "@astrojs/react";
import mdx from "@astrojs/mdx";
import svgr from "vite-plugin-svgr";

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
