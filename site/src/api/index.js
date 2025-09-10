/** whether to mock (fake) network requests */
export const mock = !!new URL(window.location.href).searchParams.has("mock");

/** general request */
export async function request(
  /** request url */
  url,
  /** fetch options */
  options = {},
  /** parse response mode */
  parse = "json",
) {
  /** make url object */
  url = new URL(url);
  /** construct request */
  const request = new Request(url, options);
  console.debug(`📞 Request ${url}`, { options, request });
  /** make request */
  const response = await fetch(request);
  /** capture error for throwing later */
  let error = "";
  /** check status code */
  if (!response.ok) error = "Response not OK";
  /** parse response */
  let parsed;
  try {
    parsed =
      parse === "text"
        ? await response.clone().text()
        : await response.clone().json();
  } catch (e) {
    error = `Couldn't parse response as ${parse}`;
  }
  console.debug(`📣 Response ${url}`, { response, parsed });
  /** throw error after details have been logged */
  if (error || parsed === undefined) throw Error(error);
  return parsed;
}
