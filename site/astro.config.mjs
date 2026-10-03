import mdx from "@astrojs/mdx";
import react from "@astrojs/react";
import sitemap from "@astrojs/sitemap";
import tailwindcss from "@tailwindcss/vite";
import { defineConfig } from "astro/config";
import svgr from "vite-plugin-svgr";
import { redirects } from "./redirects.mjs";

// https://astro.build/config
export default defineConfig({
  site: "https://strchive.org",
  integrations: [
    mdx(),
    react({}),
    sitemap({
      filter: (path) => !path.endsWith("/edit"),
    }),
    redirects(),
  ],
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
      tailwindcss(),
    ],
  },
});
