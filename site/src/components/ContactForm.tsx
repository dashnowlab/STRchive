import { createIssue } from "@/api/issue";
import Alert from "@/components/Alert";
import Button from "@/components/Button";
import Form from "@/components/Form";
import Help from "@/components/Help";
import Link from "@/components/Link";
import TextBox from "@/components/TextBox";
import { repo } from "@/layouts/meta";
import { userAgent } from "@/util/browser";
import { useQuery } from "@/util/hooks";
import { shortenUrl } from "@/util/string";
import { useLocalStorage } from "@reactuses/core";
import { IconSend } from "@tabler/icons-react";
import clsx from "clsx";
import { mapValues, startCase, truncate } from "lodash-es";

/** shared schema info for user contact info */
export const contactSchema = {
  name: {
    title: "Name",
    description: "Optional. So we know who you are.",
    placeholder: "Your Name",
  },
  username: {
    title: "GitHub",
    description: "Optional. So we can tag you and you can follow activity.",
    placeholder: "@username",
  },
  email: {
    title: "Email",
    description: "Optional. So we can contact you directly if needed.",
    placeholder: "you@example.com",
  },
};

/** contact field component props */
const contactFields = mapValues(
  contactSchema,
  ({ title, description, placeholder }) => ({
    label: (
      <>
        {title}
        <Help>{description}</Help>
      </>
    ),
    placeholder,
  }),
);

export default function ContactForm() {
  /** form state */
  let [name, setName] = useLocalStorage("contact-name", "");
  let [username, setUsername] = useLocalStorage("contact-username", "");
  let [email, setEmail] = useLocalStorage("contact-email", "");
  let [subject, setSubject] = useLocalStorage("contact-subject", "");
  let [message, setMessage] = useLocalStorage("contact-message", "");

  /** set fallbacks */
  name ||= "";
  username ||= "";
  email ||= "";
  subject ||= "";
  message ||= "";

  /** validate username */
  if (username && username.length > 0)
    username = username.replaceAll(/^@*/g, "@");

  /** extra details to include in report */
  const { pathname, search, hash } = window.location;
  const { browser, engine, os, device, cpu } = userAgent;
  const details = mapValues(
    {
      Page: [pathname + search + hash],
      Browser: [browser.name, browser.version],
      Engine: [engine.name, engine.version],
      OS: [os.name, os.version],
      Device: [device.type, device.model, device.vendor],
      CPU: [cpu.architecture],
    },
    (value) => value.filter(Boolean).join(" "),
  );

  /** message title */
  const title = truncate(
    [name.trim() || username.trim(), subject.trim()]
      .filter(Boolean)
      .join(" - "),
    { length: 250 },
  );

  /** message body */
  const body = [{ name, username, email }, details, { message }]
    .map((group) =>
      Object.entries(group)
        .map(([key, value]) => [
          `**${startCase(key)}**`,
          value.trim() ? value.trim() : "\\-",
        ])
        .flat()
        .join("\n"),
    )
    .join("\n\n");

  /** submission query */
  const {
    query: submit,
    data: response,
    status,
  } = useQuery(() =>
    createIssue({
      owner: "dashnowlab",
      repo: "STRchive",
      title,
      body,
      labels: ["contact"],
    }),
  );

  return (
    <Form onSubmit={submit}>
      <Alert type="info" className="w-0 min-w-full">
        Want to suggest a new locus or an edit to an existing one? Use the{" "}
        <Link to="/loci/new">new locus form</Link> or go to an{" "}
        <Link to="/loci#loci">existing locus</Link> and use the{" "}
        <em>suggest edit</em> form.
      </Alert>

      <div
        className={clsx(
          "grid grid-cols-[auto_minmax(--spacing(50),--spacing(150))] gap-4 max-sm:grid-cols-1 [&>label]:contents",
        )}
      >
        <TextBox {...contactFields.name} value={name} onChange={setName} />
        <TextBox
          {...contactFields.username}
          value={username}
          onChange={setUsername}
        />
        <TextBox {...contactFields.email} value={email} onChange={setEmail} />

        <TextBox
          label="Subject"
          placeholder="Subject"
          value={subject}
          onChange={setSubject}
          required
        />
        <TextBox
          label="Message"
          placeholder="Comments, questions, issues, etc."
          value={message}
          onChange={setMessage}
          multi
          required
        />
      </div>

      <Alert type={status || "info"} className="w-0 min-w-full">
        {status === "" && (
          <>
            This will make a <strong>public</strong> post on{" "}
            <Link to={repo}>our GitHub</Link> with <em>all of the above</em> and
            some debug info. You'll get a link to it once it's created.
          </>
        )}
        {startCase(status)}{" "}
        {response && (
          <Link to={response.link}>{shortenUrl(response.link)}</Link>
        )}
      </Alert>

      {status === "" && (
        <Button className="w-60 self-center" design="plain" type="submit">
          <IconSend />
          Submit
        </Button>
      )}
    </Form>
  );
}
